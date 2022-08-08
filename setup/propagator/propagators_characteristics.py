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
from matplotlib import pyplot as plt
import glob

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy import util

# Custom imports
from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust

fig1 = plt.figure(figsize=(9,4))
fig2 = plt.figure(figsize=(9,3.5))
fig3 = plt.figure(figsize=(9,3.5))
fig4 = plt.figure(figsize=(9,3.5))
base_states = None
base_dep_vars = None
base_label = None
propagators = [
    propagation_setup.propagator.cowell,
    propagation_setup.propagator.encke,
    propagation_setup.propagator.gauss_keplerian,
    propagation_setup.propagator.gauss_modified_equinoctial,
    propagation_setup.propagator.unified_state_model_quaternions,
    propagation_setup.propagator.unified_state_model_modified_rodrigues_parameters,
    propagation_setup.propagator.unified_state_model_exponential_map

]
propagators_labels = [
    "Cowell",
    "Encke",
    "Gauss-Kepler",
    "GME",
    "USM Q",
    "USM MRP",
    "USM EM"
]

for i_prop, propagator in enumerate(propagators):
    # Define problem
    SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
    SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)

    SRM_stage_2 = spherical_SRM(R_o=0.165, R_i=0.0915)
    SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)

    mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
    mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p

    # Max deflection of 5 deg
    body_fixed_thrust_direction_y = [
        [0, 0.05, 0.1, 0, 0.05],
        0   # second stage has no TVC
    ]
    body_fixed_thrust_direction_z = [
        [0, -0.05, 0.0, 0.05, 0.05],
        0   # second stage has no TVC
    ]

    MAV_ascent = ascent_framework.MAV_ascent(
        launch_epoch = 0,
        launch_lat = np.deg2rad(18.85),
        launch_lon = np.deg2rad(77.52),
        launch_h = -2.5e3,
        mass_stages = [mass_1, mass_2],
        launch_angles = [np.deg2rad(57.5), np.deg2rad(90)],
        thrust_models = [SRM_thrust_model_1, SRM_thrust_model_2],
        target_orbit_h = 300e3,
        target_orbit_i = np.deg2rad(25),
        max_a = 15 * 9.80665,
        max_AoA = np.deg2rad(4),
        body_fixed_thrust_direction_y=body_fixed_thrust_direction_y,
        body_fixed_thrust_direction_z=body_fixed_thrust_direction_z
    )

    # Setup and run simulation for both stages
    stage_res = []
    for stage in [1, 2]:
        MAV_ascent.create_bodies(
            stage=stage,
            add_sun=True,
            use_new_coeffs=True
        )
        MAV_ascent.create_accelerations(use_cpp=True, better_precision=False, thrust_fname=glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0])
        guidance_object = ascent_framework.FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
        MAV_ascent.create_initial_state()
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.keplerian_state(MAV_ascent.current_name, "Mars"))
        MAV_ascent.create_termination_settings(end_time=160*60, exact_time=True)
        MAV_ascent.create_propagator_settings(propagator_type=propagator)
        MAV_ascent.create_integrator_settings(better_accuracy=False)
        if stage == 1:
            print("Running with propagator:", propagators_labels[i_prop])
        times, states, dep_vars = MAV_ascent.run_simulation()
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,0])
        if stage == 1:
            t_b_1 = MAV_ascent.thrust.burn_time
            t_sep = times[-1]-MAV_ascent.launch_epoch
            idx_sep = len(times)-1
            if final_h < 0:
                break
        else:
            t_b_2 = MAV_ascent.thrust.burn_time

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

    label = propagators_labels[i_prop]

    states_dict = {}
    dep_vars_dict = {}
    init_state = None
    for t, state, dep_var in zip(times, states, dep_vars):
        if init_state is None:
            init_state = state
        elif np.fabs(t-t_sep) > 60:
            states_dict[t] = [np.linalg.norm(state[:3]-init_state[:3])]
            dep_vars_dict[t] = dep_var
    if base_states is None:
        base_states = states_dict
        base_dep_vars = dep_vars_dict
        base_label = propagators_labels[i_prop]
    else:
        t0, tend = 60, 157*60
        diff_times = np.arange(t0, tend, 0.1)
        dep_vars_diff = util.compare_results(base_dep_vars, dep_vars_dict, diff_times)
        dep_vars_diff_array = util.result2array(dep_vars_diff)
        dep_vars = dep_vars_diff_array[:,1:]
        states_diff = util.compare_results(base_states, states_dict, diff_times)
        states_diff_array = util.result2array(states_diff)
        states = states_diff_array[:,1:]
        label = "%s (vs %s)" % (propagators_labels[i_prop], base_label)
        eccentricities = dep_vars[:,1]
        inclinations = np.rad2deg(dep_vars[:,2])
        # positions = np.linalg.norm(states[:,:3]-states[0,:3], axis=1)
        positions = np.fabs(states[:,0])
        # angu_mom = np.cross(states[:,:3], states[:,3:6])
        # angu_mom_norm = np.linalg.norm(angu_mom, axis=1)
        # diff_times = times
        print("Final t=%.5f min, e=%.5f, i=%.5f deg, pos=%.5e m" % (times[-1]/60, eccentricities[-1], inclinations[-1], positions[-1]))
        plt.figure(fig1)
        plt.plot(diff_times/60, eccentricities, label=label, linestyle="dashed")
        plt.figure(fig2)
        plt.plot(diff_times/60, inclinations, label=label, linestyle="dashed")
        plt.figure(fig3)
        plt.plot(diff_times/60, positions/1e3, label=label, linestyle="dashed")
        # plt.figure(fig4)
        # plt.plot(diff_times/60, angu_mom_norm, label=label, linestyle="dashed")


plt.figure(fig1)
plt.xlabel("Time [min]")
plt.ylabel("Eccentricity [-]")
plt.grid()
plt.yscale("log")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/propagator_e.pdf")

plt.figure(fig2)
plt.xlabel("Time [min]")
plt.ylabel("Inclination [deg]")
plt.grid()
plt.yscale("log")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/propagator_i.pdf")

plt.figure(fig3)
plt.xlabel("Time [min]")
plt.ylabel("Error in position [m]")
plt.grid()
plt.yscale("log")
# plt.ylim(1e5, 1e7)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/propagator_pos.pdf")

plt.figure(fig4)
# mu_Mars = 4.282837e13
# circ_angu_mom = np.sqrt(mu_Mars*np.linalg.norm(states[:,:3], axis=1))
# xlims = plt.xlim()
# plt.hlines(circ_angu_mom, -1e3, 1e3, colors="black", linestyles="dashed")
# plt.xlim(xlims)
plt.xlabel("Time [min]")
plt.ylabel("Angular momentum [m/s$^2$]")
plt.grid()
plt.yscale("log")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/propagator_momentum.pdf")