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
import os

# Tudatpy imports
from tudatpy import util

# Get list of timesteps for which simulations have been run
filenames = os.listdir(sys.path[0]+"/setup/integrator/benchmark_sim_results")
filenames.remove(".gitkeep")
list_dts = sorted([float(name.replace("dt_", "").replace(".npz", "")) for name in filenames])

# print(np.diff(np.log2(list_dts))), input()

# Define function that returns the times, states, dependent variables, and number of function evaluations
# corresponding to a given integrator time step
def get_sim_results(dt):
    print("Loading state history with dt = %.4e s..." % dt)
    results = np.load(sys.path[0]+"/setup/integrator/benchmark_sim_results/dt_%.4e.npz" % dt)
    return results["times"], results["states"], results["dep_vars"], results["f_evals"]

# Get the results from the integration with the lowest timestep, to use as a baseline
baseline_dt = list_dts.pop(0)   # (remove this first dt from the list)
baseline_times, baseline_states, baseline_dep_vars, baseline_f_evals = get_sim_results(baseline_dt)
baseline_states_dict = {baseline_times[i]: baseline_states[i] for i in range(len(baseline_times))}

saved_dt = []

position_errors, velocity_errors = [], []
mass_errors = []
for i, dt in enumerate(list_dts):
    times, states, dep_vars, f_evals = get_sim_results(dt)
    states_dict = {times[i]: states[i] for i in range(len(times))}
    # Break the loop if there are less than 8 data points (not possible for Lagrange interpolator and bad results anyways)
    if len(times) < 8:
        print("Stopping, not enough data points in current propagation...")
        break
    # Compute the state difference between the baseline and the current results
    # Skip the first 10 seconds and finish one minute before to avoid interpolation artifacts
    print("Computing difference...")
    states_diff = util.compare_results(states_dict, baseline_states_dict, np.linspace(times[0]+10, min(times[0]+99*60, times[-1]), 5000))
    states_diff_array = util.result2array(states_diff)
    diff_times = states_diff_array[:,0] - states_diff_array[0,0]
    diff_pos = np.linalg.norm(states_diff_array[:,1:4], axis=1)
    diff_vel = np.linalg.norm(states_diff_array[:,4:7], axis=1)
    diff_mass = states_diff_array[:,7]
    position_errors.append(max(diff_pos)), velocity_errors.append(max(diff_vel))
    mass_errors.append(max(diff_mass))
    saved_dt.append(dt)
    # print("With a time step of %.4e [s], maximum errors: position of %.4e [m] / velocity of %.4e [m/s]" \
    #     % (dt, max(diff_pos), max(diff_vel)))
    print("With a time step of %.4e [s], final errors: position of %.4e [m] / velocity of %.4e [m/s] / mass of %.4e [kg]" \
        % (dt, position_errors[-1], velocity_errors[-1], mass_errors[-1]))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(diff_times/60, diff_pos, label="Position [m]")
    ax.plot(diff_times/60, diff_vel, label="Velocity [m/s]")
    # ax.plot(diff_mass/60, diff_vel, label="Mass [kg]")
    ax.set_xlabel("Time [min]"), ax.set_ylabel("Error")
    plt.suptitle("State error for $\Delta$t = %.3e [s] (w.r.t. $\\frac{\Delta t}{2}$)" % dt)
    ax.set_yscale("log")
    plt.grid(), plt.legend(), plt.tight_layout()
    plt.savefig(sys.path[0] + "/plots/setup/integrator/benchmark/error_dt_%.4e.pdf" % dt)
    plt.close()

    baseline_states_dict = states_dict

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(saved_dt, position_errors, marker="o", linewidth=1.5, label="Position [m]")
ax.plot(saved_dt, velocity_errors, marker="o", linewidth=1.5, label="Velocity [m/s]")
ax.plot(saved_dt, mass_errors, marker="o", linewidth=1.5, label="Mass [kg]")
ax.set_xlabel("Time step [s]"), ax.set_ylabel("Maximum error")
# ax.set_xlabel("Time step [s]"), ax.set_ylabel("Final error")
ax.set_xscale("log"), ax.set_yscale("log")
plt.grid(), plt.legend(), plt.tight_layout()
plt.savefig(sys.path[0] + "/plots/setup/integrator/benchmark/dt_vs_error.pdf")
plt.show()