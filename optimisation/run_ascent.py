import sys, os
# Add tudatpy path
if "cala" in os.getcwd():
    sys.path.append("/cala/jeremie/tudat-bundle/build/tudatpy")
else:
    sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Standard imports
import numpy as np
import pygmo as pg

# Tudat imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice
from tudatpy import util

# Custom imports
from setup.ascent_framework import MAV_ascent, FakeAeroGuidance
from thrust.models.multi_fin import multi_fin_SRM
from optimisation import ascent_problem as AP

def score_altitude(h, h_min=300e3, h_max=375e3):
    h_buffer = h_min + (h_max-h_min)/5
    if h_buffer <= h and h <= h_max:
        return 0
    elif h > h_max:
        return (h-h_max)/(450e3-h_max)
    elif h >= h_min and h < h_buffer:
        return (h_buffer-h)/(450e3-h_max)
    else:
        return (h-h_buffer)/(h_min-h_buffer)*0.75 - 9/16

def score_mass(m, m_limit=400):
    if m <= m_limit:
        return m / m_limit
    else:
        return (5*m_limit-4*m-500)/(m_limit-500)

def MAV_ascent_sim(
        thrust_angle_1,
        thrust_angle_2,
        TVC_angles_y,
        TVC_angles_z,
        thrust_model_1,
        thrust_model_2,
        return_sim_results=False
    ):
    print("Simulating MAV ascent with inputs hash", hash(str(thrust_angle_1)+str(thrust_angle_2)+str(TVC_angles_y)+str(TVC_angles_z)+str(thrust_model_1)+str(thrust_model_2)))

    mass_2 = 47.5 + thrust_model_2.M_innert + thrust_model_2.M_p
    mass_1 = 65 + mass_2 + thrust_model_1.M_innert + thrust_model_1.M_p

    ascent_model = MAV_ascent(
        launch_epoch = 0,
        launch_lat = np.deg2rad(18.85),
        launch_lon = np.deg2rad(77.52),
        launch_h = -2.5e3,
        mass_stages = [mass_1, mass_2],            
        launch_angles = [thrust_angle_1, thrust_angle_2],
        thrust_models = [thrust_model_1, thrust_model_2],
        target_orbit_h = 300e3,
        target_orbit_i = np.deg2rad(25),
        max_a = 15 * 9.80665,
        max_AoA = np.deg2rad(4),
        body_fixed_thrust_direction_y=[ TVC_angles_y, 0],
        body_fixed_thrust_direction_z=[ TVC_angles_z, 0 ]
    )

    # Setup and run simulation for both stages
    stage_res = []
    t_b_1, t_b_2 = 0, 0
    for stage in [1, 2]:
        ascent_model.create_bodies(stage=stage, add_sun=True, use_new_coeffs=True)
        ascent_model.create_accelerations(use_cpp=(stage==1), better_precision=True)
        guidance_object = FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, ascent_model.current_body, silence_warnings=True)
        ascent_model.create_initial_state()
        ascent_model.create_dependent_variables_to_save(False)
        ascent_model.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(ascent_model.current_name, "Mars"))
        ascent_model.create_termination_settings(end_time=160*60, cpu_time_termination=30)
        ascent_model.create_propagator_settings()
        ascent_model.create_integrator_settings()
        with util.redirect_std():
            times, states, dep_vars = ascent_model.run_simulation()
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,0])
        if stage == 1:
            t_b_1 = ascent_model.thrust.burn_time
            if final_h < 0:
                break
        else:
            t_b_2 = ascent_model.thrust.burn_time
    if stage == 1:
        # Extract results from first propagation only if stage 2 was not used
        times = stage_res[0][0]
        states = stage_res[0][1]
        dep_vars = stage_res[0][2]
    else:
        # Combine results from both propagations
        times = np.concatenate((stage_res[0][0], stage_res[1][0]))
        states = np.concatenate((stage_res[0][1], stage_res[1][1]))
        dep_vars = np.concatenate((stage_res[0][2], stage_res[1][2]))

    t_b = t_b_1 + t_b_2

    final_cartesian_state = states[-1,:-1]

    final_keplerian_state = element_conversion.cartesian_to_keplerian(
        final_cartesian_state,
        spice.get_body_gravitational_parameter("Mars")
    )
    final_a, final_e, final_i = final_keplerian_state[:3]
    R_Mars = spice.get_average_radius("Mars")

    h_peri, h_apo = final_a * (1 - final_e) - R_Mars, final_a * (1 + final_e) - R_Mars

    print("Final time %.2f min: apoapsis = %.2f km / periapsis = %.2f / inclination = %.2f deg with mass = %.2f kg" \
        % (times[-1]/60, h_apo/1e3, h_peri/1e3, np.rad2deg(final_i), mass_1))

    altitude_score = score_altitude(h_peri) + score_altitude(h_apo)
    mass_score = score_mass(mass_1)

    # Do some cleanup
    del(ascent_model)
    del(thrust_model_1)
    del(thrust_model_2)

    if return_sim_results:
        return altitude_score, mass_score, times, states, dep_vars

    return altitude_score, mass_score

if __name__ == "__main__":
    # Parameters
    seed = 42
    pop_size = 8
    n_generations = 5

    # Define the range in which the design variables can vary
    launch_angle_1_range = np.deg2rad([30, 60])
    launch_angle_2_range = np.deg2rad([60, 120])
    TVC_range = np.deg2rad([-5, 5])
    N_TVC_nodes = 5
    spherical_SRM_range = [[0.3, 0.2], [1.0, 0.9]]
    multi_fin_SRM_range = [[0.3, 0.1, 0.2, 0.25, 0.35, 3], [1.25, 0.285, 0.9, 0.75, 0.9, 15]]
    design_var_range = (
        [launch_angle_1_range[0], launch_angle_2_range[0], *[TVC_range[0]]*N_TVC_nodes*2, *spherical_SRM_range[0], *multi_fin_SRM_range[0]],
        [launch_angle_1_range[1], launch_angle_2_range[1], *[TVC_range[1]]*N_TVC_nodes*2, *spherical_SRM_range[1], *multi_fin_SRM_range[1]]
    )

    # Define the optimisation problem
    ascent_problem = AP.MAV_problem(design_var_range, multi_fin_SRM)
    problem = pg.problem(ascent_problem)

    # Initialise a Pygmo population
    pop = pg.population(problem, size=pop_size-1, seed=seed, b=pg.default_bfe())
    # Add a good initial guess to the population
    pop.push_back(
        x=[np.deg2rad(57.5), np.deg2rad(90), 0, 0.05, 0.1, 0, 0.05, 0, -0.05, 0.0, 0.05, 0.05, 0.6875, 0.0915/0.165, 1.05, 0.24, 0.175/0.24, 0.05/0.175, 0.02/(2*np.pi*(0.175-0.05)/20), 20],
        f=[1.7987177916981962, 0.9681385603150028])
    
    # Select the optimisation algorithm
    algo = pg.nsga2(seed=seed)
    algo.set_bfe(pg.bfe())
    algo = pg.algorithm(algo)

    print("Initial population:")
    print(pop)

    input("Press enter to start evolving the population...")

    # Run the optimisation
    for i in range(n_generations):
        print("Running generation %2d / %2d" % (i+1, n_generations))
        # Evolve the population
        pop = algo.evolve(pop)

        # Print the population
        print(pop)


# test_scoring = False
# if test_scoring:
#     import matplotlib.pyplot as plt

#     altitudes = np.arange(200e3, 500e3, 1)
#     scores = [score_altitude(h) for h in altitudes]

#     plt.plot(altitudes/1e3, scores)
#     ymin, ymax = plt.ylim()
#     plt.vlines([300, 375], -10, 10, colors=["red", "orange"])
#     plt.ylim(ymin, ymax), plt.xlim(min(altitudes)/1e3, max(altitudes)/1e3)
#     plt.xlabel("Altitude [km]"), plt.ylabel("Fitness")
#     plt.grid(), plt.tight_layout()
#     plt.show()

#     masses = np.arange(250, 550, 0.1)
#     scores = [score_mass(m) for m in masses]

#     plt.plot(masses, scores)
#     ymin, ymax = plt.ylim()
#     plt.vlines([400], -10, 10, colors=["red"])
#     plt.ylim(ymin, ymax), plt.xlim(min(masses), max(masses))
#     plt.xlabel("Mass [kg]"), plt.ylabel("Fitness")
#     plt.grid(), plt.tight_layout()
#     plt.show()