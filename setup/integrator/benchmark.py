import sys
# Add tudatpy path
sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Standard imports
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import copy

# Tudatpy imports
from tudatpy import util
from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup

# Custom imports
from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust import SRM_thrust

t0 = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17)))
SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.1725, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)
SRM_stage_2 = spherical_SRM(R_o=0.175, R_i=0.08)
SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)
mass_2 = 60+SRM_thrust_model_2.M_innert+SRM_thrust_model_2.M_p
mass_1 = 35+mass_2+SRM_thrust_model_1.M_innert+SRM_thrust_model_1.M_p
body_fixed_thrust_direction = [
    [0, 0.05, 0.1, 0, 0.05],
    0
]
# Define default ascent model
MAV_ascent_original = ascent_framework.MAV_ascent(
    launch_epoch = t0,    # MAV-­LL­-01
    launch_lat = np.deg2rad(18.85),     # MAV­-LL-­03
    launch_lon = np.deg2rad(77.52),     # MAV­-LL-­03
    launch_h = -2.5e3,                  # MAV­-LL-­04
    mass_stages = [mass_1, mass_2],            
    launch_angles = [np.deg2rad(57.5), np.deg2rad(90)],       # MAV­-LL-­06 + guesstimate # angle is w.r.t vertical
    thrust_models = [SRM_thrust_model_1, SRM_thrust_model_2],
    target_orbit_h = 300e3,             # MAV­-OSO­-01
    target_orbit_i = np.deg2rad(25),    # MAV­-OSO­-03
    max_a = 15 * 9.80665,               # MAV­-LL-­02
    max_AoA = np.deg2rad(4),            # MAV­-LL-­05
    body_fixed_thrust_direction=body_fixed_thrust_direction
)

def run_all(dt):
    MAV_ascent = copy.deepcopy(MAV_ascent_original)
    stage_res = []
    # Setup and run simulation for stage 1 then 2
    print("Running with dt = %.3e s" % dt)
    # f_evals_list = []
    for stage in [1, 2]:
        MAV_ascent.create_bodies(stage=stage)
        MAV_ascent.create_accelerations()
        guidance_object = ascent_framework.FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
        MAV_ascent.create_initial_state()
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(MAV_ascent.current_name, "Mars"))
        MAV_ascent.create_termination_settings(end_time=95*60)
        MAV_ascent.create_propagator_settings()
        MAV_ascent.create_integrator_settings(fixed_step=dt)
        with util.redirect_std():
            times, states, dep_vars, f_evals = MAV_ascent.run_simulation(return_count=True)
        stage_res.append([times, states, dep_vars, f_evals])

    # Combine and save results from both propagations
    if stage == 1:
        # Extract results from first propagation only if stage 2 was not used
        times = stage_res[0][0]
        states = stage_res[0][1]
        dep_vars = stage_res[0][2]
        f_evals = stage_res[0][3]
    else:
        # Combine results from both propagations
        times = np.concatenate((stage_res[0][0], stage_res[1][0]))
        states = np.concatenate((stage_res[0][1], stage_res[1][1]))
        dep_vars = np.concatenate((stage_res[0][2], stage_res[1][2]))
        f_evals = stage_res[0][3] + stage_res[1][3]

    np.savez("setup/integrator/benchmark_sim_results/dt_%.4e"%dt, times=times, states=states, dep_vars=dep_vars, f_evals=f_evals)
    print(" -> %.3e function evaluations." % f_evals)
    

max_pos_diffs, max_vel_diffs = [], []
dts = []
# smallest_dt, largest_dt = 5e-1*1.01, 10 # 1.01 try to avoid perfect numerical values (without decimals)

# dt = smallest_dt
start_dt = 12.45 # (do no use perfect int to avoid artifacts with perfect numerical values)
dt = start_dt
# more_accurate_states, more_accurate_altitudes = run_all(dt/2)
while True:#dt <= largest_dt*2:
    # current_states, current_altitudes = run_all(dt)
    run_all(dt)

    # print("current:", list(current_states.values())[50])
    # print("better:", list(more_accurate_states.values())[50])
    # benchmark_diff = util.compare_results(current_states, more_accurate_states, np.linspace(t0+10, t0+19*60, 5000))
    # benchmark_diff_array = util.result2array(benchmark_diff)
    # diff_times = benchmark_diff_array[:,0] - benchmark_diff_array[0,0]
    # diff_pos = np.linalg.norm(benchmark_diff_array[:,1:4], axis=1)
    # diff_vel = np.linalg.norm(benchmark_diff_array[:,4:7], axis=1)
    # # plt.plot(diff_times/60, diff_pos, label="Position [m]")
    # # plt.plot(diff_times/60, diff_vel, label="Velocity [m/s]")
    # more_accurate_states_array = util.result2array(more_accurate_states)
    # current_states_array = util.result2array(current_states)
    # # plt.plot((more_accurate_states_array[:,0]-more_accurate_states_array[0,0])/60, more_accurate_altitudes/1e3, label="dt = %.2e [s]" % (dt/2))
    # # plt.plot((current_states_array[:,0]-current_states_array[0,0])/60, current_altitudes/1e3, label="dt = %.2e [s]" % (dt))
    # # plt.xlabel("Time [min]"), plt.ylabel("Altitude [km]")
    # # plt.grid(), plt.legend(), plt.tight_layout()

    # max_pos_diffs.append(max(diff_pos)), max_vel_diffs.append(max(diff_vel))
    # print("Maximum position error of %.3e m, velocity error of %.3e m/s" % (max_pos_diffs[-1], max_vel_diffs[-1]))
    # dts.append(dt)
    # more_accurate_states = current_states
    # more_accurate_altitudes = current_altitudes
    # plt.show()
    dt /= 2

    # fig = plt.figure(figsize=(7, 6))
    # ax = fig.add_subplot(111, projection="3d")
    # ax.plot(more_accurate_states_array[:,1], more_accurate_states_array[:,2], more_accurate_states_array[:,3], label="dt = %.2e [s]" % (dt/2))
    # ax.plot(current_states_array[:,1], current_states_array[:,2], current_states_array[:,3], label="dt = %.2e [s]" % (dt))
    # ax.legend()
    # plt.show()

plt.plot(dts, max_pos_diffs, label="Position [m]", linewidth=1, markersize=4, marker="o")
plt.plot(dts, max_vel_diffs, label="Velocity [m/s]", linewidth=1, markersize=4, marker="o")
plt.xlabel("Step size [s]"), plt.ylabel("Error")
plt.yscale("log"), plt.xscale("log")
plt.grid(), plt.legend(), plt.tight_layout()
plt.show()

# # Save the benchmark values
# np.savez(sys.path[0]+"/setup/integrator/benchmark", times=times, states=states, f_evals=f_evals)