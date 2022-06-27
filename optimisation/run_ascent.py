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
import time as T
import sqlite3

# Tudat imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice
from tudatpy import util
from tudatpy.kernel.math import interpolators

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
        print_times=False,
        save_to_db=None
    ):
    # print()
    # print("Simulating MAV ascent with inputs:")#, hash(str(thrust_angle_1)+str(thrust_angle_2)+str(TVC_angles_y)+str(TVC_angles_z)+str(thrust_model_1)+str(thrust_model_2)))
    # print("Thrust model 1:", thrust_model_1.geometry_model)
    # print("Thrust model 2:", thrust_model_2.geometry_model)
    # print("Thrust angle 1:", thrust_angle_1)
    # print("Thrust angle 2:", thrust_angle_2)
    # print("TVC angle y:", TVC_angles_y)
    # print("TVC angle z:", TVC_angles_z)
    # print()

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
        ascent_model.create_bodies(stage=stage, add_sun=True, use_new_coeffs=True, custom_exponential_model=True)
        t0 = T.time()
        try:
            ascent_model.create_accelerations(use_cpp=(stage==1), better_precision=True)
        except ValueError:
            return 100, 100
        guidance_object = FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, ascent_model.current_body, silence_warnings=True)
        ascent_model.create_initial_state()
        ascent_model.create_dependent_variables_to_save(False)
        ascent_model.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(ascent_model.current_name, "Mars"))
        ascent_model.create_termination_settings(end_time=160*60, cpu_time_termination=30)
        ascent_model.create_propagator_settings()
        ascent_model.create_integrator_settings()
        with util.redirect_std():
            t1 = T.time()
            times, states, dep_vars = ascent_model.run_simulation()
            t2 = T.time()
        dt = t2-t1
        if print_times:
            print("Burn simulation took", t1-t0, "seconds")
            print("Simulation took", dt, "seconds")
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,0])
        if stage == 1:
            t_b_1 = ascent_model.thrust.burn_time
            t_sep = times[-1]
            if final_h < 0 or dt > 29:
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

    final_cartesian_state = states[-1,:-1]

    final_keplerian_state = element_conversion.cartesian_to_keplerian(
        final_cartesian_state,
        spice.get_body_gravitational_parameter("Mars")
    )
    final_a, final_e, final_i = final_keplerian_state[:3]
    R_Mars = spice.get_average_radius("Mars")

    h_peri, h_apo = final_a * (1 - final_e) - R_Mars, final_a * (1 + final_e) - R_Mars

    add_print = ""
    if save_to_db is not None:
        add_print = " (saved with id %i)" % save_to_db

    print("Final time %.2f min: apoapsis = %.2f km / periapsis = %.2f / inclination = %.2f deg with mass = %.2f kg%s" \
        % (times[-1]/60, h_apo/1e3, h_peri/1e3, np.rad2deg(final_i), mass_1, add_print))

    # Resample results
    if save_to_db is not None:
        if len(times) > 50:
            times_resampled_1 = np.linspace(0, t_sep, num=1030)[:-30]
            times_resampled_2 = np.linspace(t_sep, times[-1], num=1030)[30:]
            times_resampled = np.concatenate((times_resampled_1, times_resampled_2))
            states_dict = {epoch: state for epoch, state in zip(times, states)}
            dep_vars_dict = {epoch: dep_var for epoch, dep_var in zip(times, dep_vars)}
            interpolator_settings = interpolators.lagrange_interpolation(8, boundary_interpolation=interpolators.use_boundary_value)
            states_interpolator = interpolators.create_one_dimensional_vector_interpolator(states_dict, interpolator_settings)
            dep_vars_interpolator = interpolators.create_one_dimensional_vector_interpolator(dep_vars_dict, interpolator_settings)
            states_resampled = np.asarray([states_interpolator.interpolate(epoch) for epoch in times_resampled])
            dep_vars_resampled = np.asarray([dep_vars_interpolator.interpolate(epoch) for epoch in times_resampled])
        else:
            times_resampled = times
            states_resampled = states
            dep_vars_resampled = dep_vars

    h_p_score, h_a_score = score_altitude(h_peri), score_altitude(h_apo)
    mass_score = score_mass(mass_1)

    if save_to_db is not None:
        # Connect to the database
        con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
        cur = con.cursor()
        # Insert results
        req = "UPDATE solutions_multi_fin "
        # req += "(h_p_score, h_a_score, mass_score, final_time, h_a, h_p, mass, inclination, t_b_1, t_b_2)"
        req += "SET h_p_score = ?, h_a_score = ?, mass_score = ?, final_time = ?, h_a = ?, h_p = ?, mass = ?, inclination = ?, t_b_1 = ?, t_b_2 = ? "
        req += "WHERE id = ?"
        cur.execute(req, (h_p_score, h_a_score, mass_score, times[-1], h_apo, h_peri, mass_1, final_i, t_b_1, t_b_2, save_to_db))
        # Close connection
        con.commit()
        con.close()
        # Save resampled results to file
        np.savez(sys.path[0]+"/optimisation/sim_results/%i.npz" % save_to_db, times=times_resampled, states=states_resampled, dep_vars=dep_vars_resampled)


    # Do some cleanup
    del(ascent_model)
    del(thrust_model_1)
    del(thrust_model_2)

    return h_p_score, h_a_score, mass_score


if __name__ == "__main__":
    # Parameters
    seed = 42
    pop_size = 8*9
    n_generations = 30
    tuning = False
    param_tuning = False
    
    # Select the optimisation algorithm
    # BFE: GACO, MACO, NSGA2, NSPSO, PSO_GEN
    # Multi-objectives: IHS, NSGA2, MOEAD, MACO, NSPSO
    # Turn multi-objective to single-objective: see https://esa.github.io/pygmo2/problems.html#pygmo.decompose
    algo_name = "NSPSO" # NSGA2, MACO, NSPSO

    # NSPSO algorithm parameters
    diversity_mechanism = "niche count" # "crowding distance" or "niche count" or "max min"
    leader_selection_range = 5 # < 100
    omega   = 0.8   # in ]0,1[
    v_coeff = 0.05   # in ]0,1]
    chi     = 0.25   # > 0
    c1      = 0.01  # > 0
    c2      = 0.25   # > 0

    if param_tuning:
        param_to_change = {
            "diversity_mechanism": ["crowding distance", "niche count", "max min"],
            "leader_selection_range": [1, 2, 3, 5, 10],
            "omega": [0.4, 0.5, 0.6, 0.7, 0.8],
            "v_coeff": [0.05, 0.15, 0.25, 0.5, 0.75],
            "chi": [0.1, 0.25, 0.5, 1.5, 3.0],
            "c1": [0.01, 0.05, 0.15, 0.3, 0.5, 1.0],
            "c2": [0.1, 0.25, 0.5, 1.0]
        }
        for idx_param_to_change in range(len(param_to_change)):
            param_name, param_vals = list(param_to_change.items())[idx_param_to_change]
            for param_val in param_vals:
                if type(param_val) == str:
                    exec("%s = '%s'" % (param_name, param_val))
                else:
                    exec("%s = %s" % (param_name, param_val))
                param_tune_name = "%s_%s" % (param_name, param_val)

                msg = "Running optimisation with %s = %s" % (param_name, param_val)
                print("*"*(len(msg)+12))
                print("***** %s *****" % msg)
                print("*"*(len(msg)+12))
            
    ##################################################
    ### Indent three times to run parameter tuning ###
    ##################################################

    # Define the range in which the design variables can vary
    launch_angle_1_range = np.deg2rad([47.5, 60])
    launch_angle_2_range = np.deg2rad([70, 110])
    TVC_range = np.deg2rad([-5, 5])
    N_TVC_nodes = 5
    spherical_SRM_range = [[0.3, 0.2], [1.0, 0.9]]
    multi_fin_SRM_range = [[0.3, 0.1, 0.2, 0.25, 0.35, 3], [1.25, 0.285, 0.9, 0.75, 0.9, 20]]
    design_var_range = (
        [launch_angle_1_range[0], launch_angle_2_range[0], *[TVC_range[0]]*N_TVC_nodes, *spherical_SRM_range[0], *multi_fin_SRM_range[0]],
        [launch_angle_1_range[1], launch_angle_2_range[1], *[TVC_range[1]]*N_TVC_nodes, *spherical_SRM_range[1], *multi_fin_SRM_range[1]]
    )

    # Define the optimisation problem
    if tuning:
        ascent_problem = AP.MAV_problem(design_var_range, multi_fin_SRM, save_to_db="%s_%s_tuning_~N~" % (algo_name, seed))
    elif param_tuning:
        ascent_problem = AP.MAV_problem(design_var_range, multi_fin_SRM, save_to_db="%s_%s_pt_~N~" % (param_tune_name, seed))
    else:
        ascent_problem = AP.MAV_problem(design_var_range, multi_fin_SRM)
    problem = pg.problem(ascent_problem)

    # Initialise a Pygmo population
    pop = pg.population(problem, size=0, seed=seed, b=pg.default_bfe())

    # Set the algorithm
    if algo_name == "NSGA2":
        algo = pg.nsga2(seed=seed)
    elif algo_name == "MACO":
        algo = pg.maco(ker=pop_size*2//3, seed=seed, memory=True)
    elif algo_name == "NSPSO":
        algo = pg.nspso(seed=seed, omega=omega, c1=c1, c2=c2, chi=chi, v_coeff=v_coeff, leader_selection_range=leader_selection_range, diversity_mechanism=diversity_mechanism, memory=True)
    else:
        raise ValueError("Unknown algorithm")
    algo.set_bfe(pg.bfe())
    algo = pg.algorithm(algo)

    # Connect to database
    con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
    cur = con.cursor()

    # Set best results from design space exploration as initial population
    samples_weights = {
        "init_angle_only": 2,
        "TVC_only": 1,
        "SRM_only": 3,
        "all": 1
    }
    sum_weights = sum(samples_weights.values())
    samples_weights = {k: v/sum_weights for k, v in samples_weights.items()}
    ids = []
    for dv_used, fraction in samples_weights.items():
        sim_num = int(pop_size*fraction)
        if dv_used == "all":
            sim_num = pop_size - len(ids)
        req = "SELECT * FROM solutions_multi_fin WHERE dv_used = '%s' ORDER BY h_a_score+h_p_score+mass_score ASC LIMIT %i"%(dv_used, sim_num)
        cur.execute(req)
        res = cur.fetchall()
        [ids.append(row[0]) for row in res]
        for row in res:
            # Add to population
            pop.push_back(
                x=row[1:3]+row[8:21],
                f=row[21:24])
    con.close()

    # print("Initial population:")
    # print(pop)

    # Run the optimisation
    for i in range(1, n_generations+1):
        print("*** Running generation %2d / %2d with %s (seed %s) ***" % (i, n_generations, algo_name, seed))

        # Evolve the population
        pop = algo.evolve(pop)

        # # Print the population
        # print(pop)

        
    ##################################################
    ### Indent three times to run parameter tuning ###
    ##################################################




# test_scoring = False
# if test_scoring:
#     import matplotlib.pyplot as plt

#     altitudes = np.arange(-200e3, 2000e3, 1)
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
#     plt.vlines([400], -10, 10, colors=["red"])
#     plt.ylim(ymin, ymax), plt.xlim(min(masses), max(masses))
#     plt.xlabel("Mass [kg]"), plt.ylabel("Fitness")
#     plt.grid(), plt.tight_layout()
#     plt.show()
