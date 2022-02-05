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

# Tudatpy imports
from tudatpy import util

# Define default ascent model
import ascent
MAV_ascent = ascent.MAV_ascent(
    launch_epoch = 0,
    launch_lat = 0*np.deg2rad(18.85),     # MAV­-LL-­03
    launch_lon = 0*np.deg2rad(77.52),     # MAV­-LL-­03
    launch_h = -0*2.5e3,                  # MAV­-LL-­04
    mass_1 = [370, 185],                # MAV-­VM­-03 + LS p.10
    mass_2 = [80, 40],                  # LS p.10
    launch_angles = [np.deg2rad(60), np.deg2rad(90)], # MAV­-LL-­06 + guesstimate # angle is w.r.t vertical
    thrust_magnitudes = [9750, 6750],   # adapted from LS p.10
    target_orbit_h = 300e3,             # MAV­-OSO­-01
    target_orbit_i = np.deg2rad(25),    # MAV­-OSO­-03
    max_a = 15 * 9.80665,               # MAV­-LL-­02
    max_AoA = np.deg2rad(4),            # MAV­-LL-­05
    staging_altitude = 137.5e3
)

def run_all(dt):
    # Setup and run simulation for stage 1 then 2
    # times_list = []
    states_list = []
    print("Running with dt = %.3e s" % dt)
    # f_evals_list = []
    for stage in [1, 2]:
        MAV_ascent.create_bodies(stage=stage)
        MAV_ascent.create_accelerations()
        MAV_ascent.create_initial_state()
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.create_termination_settings(end_time=75*60)
        MAV_ascent.create_propagator_settings()
        MAV_ascent.create_integrator_settings(fixed_step=dt)
        with util.redirect_std():
            states, dep_vars = MAV_ascent.run_simulation(return_raw=True)
        states_list.append(states)
        #times, states, dep_vars, f_evals = MAV_ascent.run_simulation(return_count=True)
        #times_list.append(times)
        #states_list.append(states)
        #f_evals_list.append(f_evals)

    # Combine and return results from both propagations
    return {**states_list[0], **states_list[1]}
    # times = np.concatenate((times_list[0], times_list[1]))
    # states = np.concatenate((states_list[0], states_list[1]))

max_pos_diffs, max_vel_diffs = [], []
dts = []
smallest_dt, largest_dt = 1e-4*1.01, 1 # 1.01 try to avoid perfect numerical values (without decimals)
dt = smallest_dt
more_accurate_states = run_all(dt/2)
while dt <= largest_dt*2:
    current_states = run_all(dt)
    benchmark_diff = util.compare_results(current_states, more_accurate_states, np.linspace(10, 75*60-10, 5000))
    benchmark_diff_array = util.result2array(benchmark_diff)
    diff_times = benchmark_diff_array[:,0]
    diff_pos = np.linalg.norm(benchmark_diff_array[:,1:4], axis=1)
    diff_vel = np.linalg.norm(benchmark_diff_array[:,4:7], axis=1)
    # plt.plot(diff_times/60, diff_pos, label="Position [m]")
    # plt.plot(diff_times/60, diff_vel, label="Velocity [m/s]")
    # plt.xlabel("Time [min]"), plt.ylabel("Error")
    # plt.yscale("log")
    # plt.grid(), plt.legend(), plt.tight_layout()
    # plt.show()
    max_pos_diffs.append(max(diff_pos)), max_vel_diffs.append(max(diff_vel))
    print("Maximum position error of %.3e m, velocity error of %.3e m/s" % (max_pos_diffs[-1], max_vel_diffs[-1]))
    dts.append(dt)
    more_accurate_states = current_states
    dt *= 2

plt.plot(dts, max_pos_diffs, label="Position [m]", linewidth=1, markersize=4, marker="o")
plt.plot(dts, max_vel_diffs, label="Velocity [m/s]", linewidth=1, markersize=4, marker="o")
plt.xlabel("Step size [s]"), plt.ylabel("Error")
plt.yscale("log"), plt.xscale("log")
plt.grid(), plt.legend(), plt.tight_layout()
plt.show()

# # Save the benchmark values
# np.savez(sys.path[0]+"/setup/integrator/benchmark", times=times, states=states, f_evals=f_evals)