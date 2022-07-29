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
from datetime import datetime

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import time_conversion
from tudatpy import util

# Custom imports
from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust

# Parameters
analyse_gravity = False
analyse_atmo = True
analyse_misc = False
analyse_other_pm = False

n_accs = 0
if analyse_gravity:
    n_accs = 7
    f_name = "gravity"
elif analyse_atmo:
    n_accs = 3
    f_name = "atmo"
elif analyse_misc:
    n_accs = 3
    f_name = "misc"
elif analyse_other_pm:
    n_accs = 6
    f_name = "other_pm"
allowable_errors = [5e3-3.6e3, 5-3.3]

def environment_names_and_settings(i, acc_dict):
    if analyse_gravity:
        acc_dict["Mars"] = [propagation_setup.acceleration.aerodynamic()]
        if i == 0:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(14, 14))
            name = "SH D/O 14"
        elif i == 1:
            acc_dict["Mars"].append(propagation_setup.acceleration.point_mass_gravity())
            name = "Point Mass"
        elif i == 2:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(2, 2))
            name = "SH D/O 2"
        elif i == 3:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(4, 4))
            name = "SH D/O 4"
        elif i == 4:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(6, 6))
            name = "SH D/O 6"
        elif i == 5:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(8, 8))
            name = "SH D/O 8"
        elif i == 6:
            acc_dict["Mars"].append(propagation_setup.acceleration.spherical_harmonic_gravity(10, 10))
            name = "SH D/O 10"
    elif analyse_atmo:
        acc_dict["Mars"] = [
            propagation_setup.acceleration.aerodynamic(),
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 4)
        ]
        if i == 0:
            name = "Exponential atmosphere"
        elif i == 1:
            name = "Two-steps exponential atmosphere"
        # elif i == 2:
        #     name = "GRAM 2010"
        elif i == 2:
            name = "MCD 5.3"
    elif analyse_misc:
        acc_dict["Mars"] = [
            propagation_setup.acceleration.aerodynamic(),
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 4)
        ]
        if i == 0:
            name = "None"
        elif i == 1:
            acc_dict["Mars"].append(propagation_setup.acceleration.relativistic_correction(use_schwarzschild = True))
            name = "Schwarzschild relativistic correction"
        elif i == 2:
            acc_dict["Sun"] = [propagation_setup.acceleration.cannonball_radiation_pressure()]
            name = "Radiation pressure"
    elif analyse_other_pm:
        bodies = ["Mars", "Sun", "Jupiter", "Venus", "Earth", "Saturn"]
        acc_dict["Mars"] = [
            propagation_setup.acceleration.aerodynamic(),
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 4)
        ]
        if i == 0:
            name = "None"
        else:
            acc_dict[bodies[i]] = [propagation_setup.acceleration.point_mass_gravity()]
            name = "%s PM" % bodies[i]
    return name, acc_dict

max_h = None
fig1 = plt.figure(figsize=(9,5))
fig2 = plt.figure(figsize=(9,5))
if analyse_atmo:
    fig3 = plt.figure(figsize=(9,3.5))
base_states = None
base_label = None
for i_acc in range(n_accs):
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
    epoch = 0#time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 1, 1)))
    MAV_ascent = ascent_framework.MAV_ascent(
        launch_epoch = epoch if analyse_other_pm else 0,
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
            include_radiation_pressure=analyse_misc,
            more_bodies=analyse_other_pm,
            custom_exponential_model=(analyse_atmo and i_acc==1),
            use_MCD=(analyse_atmo and i_acc==2),
            # use_GRAM=(analyse_atmo and i_acc==2)
        )
        accelerations_dict = MAV_ascent.create_accelerations(use_cpp=(stage==1), only_thrust_dict=True, thrust_fname=glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0])
        label_acc, accelerations_dict = environment_names_and_settings(i_acc, accelerations_dict)
        MAV_ascent.acceleration_models = propagation_setup.create_acceleration_models(
            MAV_ascent.bodies,
            {MAV_ascent.current_name: accelerations_dict},
            MAV_ascent.bodies_to_propagate,
            MAV_ascent.central_bodies
        )
        if stage == 1:
            print("Running with environment accelerations %i/%i: %s..." % (i_acc+1, n_accs, label_acc))
        guidance_object = ascent_framework.FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
        MAV_ascent.create_initial_state()
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(MAV_ascent.current_name, "Mars"))
        MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.density(MAV_ascent.current_name, "Mars"))
        MAV_ascent.create_termination_settings(end_time=160*60, exact_time=True)
        MAV_ascent.create_propagator_settings()
        MAV_ascent.create_integrator_settings()
        # print("Running the simulation for stage %i" % stage)
        times, states, dep_vars = MAV_ascent.run_simulation()
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,0])
        # print("Max altitude at end of stage %i at %.2f min of %.2f km, vehicle mass of %.2f kg..." % \
        #     (stage, (times[-1]-MAV_ascent.launch_epoch)/60, final_h/1e3, states[-1,-1]))
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

    # Plot position over time
    states_dict = {}
    for t, state in zip(times, states):
        if np.fabs(t-t_sep) > 60:
            states_dict[t] = state
    if base_states is None and not analyse_atmo:
        base_states = states_dict
        base_label = label_acc
    else:
        if not analyse_atmo:
            t0, tend = 60, 157*60
            states_diff_1 = util.compare_results(base_states, states_dict, np.arange(t0, t_sep-90, 0.1))
            states_diff_2 = util.compare_results(base_states, states_dict, np.arange(t_sep+90, tend, 0.1))
            states_diff_array_1 = util.result2array(states_diff_1)
            states_diff_array_2 = util.result2array(states_diff_2)
            states_diff_array = np.concatenate((states_diff_array_1, states_diff_array_2))
            diff_times = states_diff_array[:,0]
            positions_diff = np.linalg.norm(states_diff_array[:,1:4], axis=1)
            velocity_diff = np.linalg.norm(states_diff_array[:,4:7], axis=1)
            label = label_acc if base_label == "None" else "%s (vs %s)" % (label_acc, base_label)
        else:
            diff_times = times
            positions_diff = np.linalg.norm(states[:,:3]-states[0,:3], axis=1)
            velocity_diff = np.linalg.norm(states[:,3:]-states[0,3:], axis=1)
            label = label_acc
        plt.figure(fig1)
        plt.plot(diff_times/60, positions_diff, label=label)
        plt.figure(fig2)
        plt.plot(diff_times/60, velocity_diff, label=label)
        print("For %s, final difference of %.2e m and %.2e m/s" %(label, positions_diff[-1], velocity_diff[-1]))
    if analyse_atmo:
        plt.figure(fig3)
        idx_crop = np.argmax(dep_vars[:,0])
        # Plot altitude vs density
        plt.plot(dep_vars[:idx_crop,0]/1e3, dep_vars[:idx_crop,1], label=label_acc)
        if i_acc == 2:
            max_h = max(dep_vars[:,0])/1e3


plt.figure(fig1)
if not analyse_atmo:
    xlims = plt.xlim()
    plt.hlines(allowable_errors[0], -1e3, 1e3, colors="orange", linestyles="dashed")
    plt.xlim(xlims)
else:
    plt.ylim(1e4, 3e7)
plt.legend(loc="lower right")
plt.xlabel("Time [min]")
plt.ylabel("$r(t) - r(0)$ [m]")
plt.title("Position over time")
plt.grid()
plt.yscale("log")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/environment/accelerations_%s_position.pdf" % f_name)

plt.figure(fig2)
if not analyse_atmo:
    xlims = plt.xlim()
    plt.hlines(allowable_errors[1], -1e3, 1e3, colors="orange", linestyles="dashed")
    plt.xlim(xlims)
else:
    plt.ylim(9e2, 1e4)
plt.xlabel("Time [min]")
plt.ylabel("$v(t) - v(0)$ [m/s]")
plt.title("Velocity over time")
plt.grid()
plt.yscale("log")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/environment/accelerations_%s_velocity.pdf" % f_name)

if analyse_atmo:
    plt.figure(fig3)
    plt.xlabel("Altitude [km]")
    plt.ylabel("Density [kg/m3]")
    plt.xlim(-5, max_h)
    plt.grid()
    plt.yscale("log")
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/setup/environment/accelerations_%s_density.pdf" % f_name)