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
plt.switch_backend('agg')
import os
import sqlite3

# Tudatpy imports
from tudatpy import util

current_stage = 1
powered = True
only_thrust = False
thrust_rk4 = True

if only_thrust:
    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/setup/integrator/integrator_tuning/database.db")
    cur = con.cursor()

# Get list of timesteps for which simulations have been run
filenames = os.listdir(sys.path[0]+"/setup/integrator/benchmark_sim_results")
list_dts = []
for name in filenames:
    try:
        if only_thrust:
            list_dts.append(float(name.replace("thrust_%s_%i_dt_"%("rk4" if thrust_rk4 else "euler", current_stage), "").replace(".npz", "")))
        else:
            list_dts.append(float(name.replace("%i_%s_dt_"%(current_stage, "V" if powered else "X"), "").replace(".npz", "")))
    except ValueError:
        pass
list_dts = sorted(list_dts)

# Define function that returns the times, states, dependent variables, and number of function evaluations
# corresponding to a given integrator time step
def get_sim_results(dt):
    if only_thrust:
        results = np.load(sys.path[0]+"/setup/integrator/benchmark_sim_results/thrust_%s_%i_dt_%.4e.npz"%("rk4" if thrust_rk4 else "euler", current_stage, dt))
        try:
            return results["times"], results["masses"], results["magnitudes"], results["fevals"]
        except KeyError:
            return results["times"], results["masses"], results["magnitudes"], 0
    else:
        results = np.load(sys.path[0]+"/setup/integrator/benchmark_sim_results/%i_%s_dt_%.4e.npz"%(current_stage, "V" if powered else "X", dt))
        return results["times"], results["states"], results["dep_vars"], results["f_evals"]

saved_dt = []

position_errors, velocity_errors = [], []
mass_errors, magnitude_errors = [], []
fevals = []
for i, dt in enumerate(list_dts):
    if dt > 5:
        break
    # Find the first time step twice as low as the current one
    try:
        baseline_dt_idx = np.where(np.asarray(list_dts) < dt/2)[0][-1]
        baseline_dt = list_dts[baseline_dt_idx]
    except IndexError:
        baseline_dt = None
    #baseline_dt = 3.1465e-8

    if baseline_dt is not None and baseline_dt != dt:
        if only_thrust:
            times, masses, magnitudes, f_evals = get_sim_results(dt)
            states_dict = {times[i]: [masses[i], magnitudes[i]] for i in range(len(times))}

            baseline_times, baseline_masses, baseline_magnitudes, _ = get_sim_results(baseline_dt)
            baseline_states_dict = {baseline_times[i]: [baseline_masses[i], baseline_magnitudes[i]] for i in range(len(baseline_times))}
        else:
            try:
                times, states, dep_vars, f_evals = get_sim_results(dt)

                baseline_times, baseline_states, baseline_dep_vars, baseline_f_evals = get_sim_results(baseline_dt)
            except ValueError:
                print("Stopping, not enough data points in current propagation...")
                break
            if times is None or dt > 0.78:
                print("Stopping, not enough data points in current propagation...")
                break
            states_dict = {times[i]: states[i] for i in range(len(times))}

            baseline_states_dict = {baseline_times[i]: baseline_states[i] for i in range(len(baseline_times))}
        # Break the loop if there are less than 8 data points (not possible for Lagrange interpolator and bad results anyways)
        if len(times) < 8:
            print("Stopping, not enough data points in current propagation...")
            break
        # Compute the state difference between the baseline and the current results
        # Skip the first and last 5 steps to avoid interpolation artifacts
        min_time = max(times[5]+5*dt, baseline_times[5]+5*baseline_dt)
        max_time = min(times[-5]-5*dt, baseline_times[-5]-5*baseline_dt)
        states_diff = util.compare_results(states_dict, baseline_states_dict, np.linspace(min_time, max_time, 500))
        states_diff_array = util.result2array(states_diff)
        diff_times = states_diff_array[:,0] - states_diff_array[0,0]
        if only_thrust:
            diff_mass = np.fabs(states_diff_array[:,1])
            diff_mag = np.fabs(states_diff_array[:,2])
            magnitude_errors.append(diff_mag[-1])
        else:
            diff_pos = np.linalg.norm(np.fabs(states_diff_array[:,1:4]), axis=1)
            diff_vel = np.linalg.norm(np.fabs(states_diff_array[:,4:7]), axis=1)
            if powered:
                diff_mass = np.fabs(states_diff_array[:,7])
            position_errors.append(diff_pos[-1]), velocity_errors.append(diff_vel[-1])
            # position_errors.append(max(diff_pos)), velocity_errors.append(max(diff_vel))
        if powered:
            mass_errors.append(diff_mass[-1])
            # mass_errors.append(max(diff_mass))
        saved_dt.append(dt)
        fevals.append(f_evals)
        if only_thrust:
            print("With a time step of %.4e [s], stage %i, final errors: thrust magnitude of %.4e [N] / mass of %.4e [kg] with %.4e f evals" \
                % (dt, current_stage, magnitude_errors[-1], mass_errors[-1], float(f_evals)))
            res = cur.execute("SELECT * FROM integrator_thrust WHERE stage=? AND method=? AND dt=?", (current_stage, "rk4" if thrust_rk4 else "euler", dt))
            if res.fetchone() is None:
                cur.execute("INSERT INTO integrator_thrust (method, f_evals, error_thrust, error_mass, stage, dt) VALUES (?, ?, ?, ?, ?, ?)", \
                    ("rk4" if thrust_rk4 else "euler", float(f_evals), magnitude_errors[-1], mass_errors[-1], current_stage, dt))
            else:
                cur.execute("UPDATE integrator_thrust SET f_evals=?, error_thrust=?, error_mass=? WHERE method=? AND dt=? AND stage=?", \
                    (float(f_evals), magnitude_errors[-1], mass_errors[-1], "rk4" if thrust_rk4 else "euler", dt, current_stage))
            con.commit()
        elif powered:
            print("With a time step of %.4e [s], stage %i, %s, final errors: position of %.4e [m] / velocity of %.4e [m/s] / mass of %.4e [kg] with %.4e f evals" \
                % (dt, current_stage, "powered" if powered else "unpowered", position_errors[-1], velocity_errors[-1], mass_errors[-1], f_evals))
        else:
            print("With a time step of %.4e [s], stage %i, %s, final errors: position of %.4e [m] / velocity of %.4e [m/s] with %.4e f evals" \
                % (dt, current_stage, "powered" if powered else "unpowered", position_errors[-1], velocity_errors[-1], f_evals))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        if only_thrust:
            ax.plot(diff_times, diff_mag, label="Magnitude [N]")
            ax.plot(diff_times, diff_mass, label="Mass [kg]")
            ax.set_xlabel("Time [s]")
        else:
            ax.plot(diff_times/60, diff_pos, label="Position [m]")
            ax.plot(diff_times/60, diff_vel, label="Velocity [m/s]")
            if powered:
                ax.plot(diff_times/60, diff_mass, label="Mass [kg]")
            ax.set_xlabel("Time [min]")
        ax.set_ylabel("Final error")
        ax.set_yscale("log")
        plt.suptitle("State error for $\Delta$t = %.3e [s] (w.r.t. $\Delta$t = %.3e [s], %s)" % (dt, baseline_dt, "RK 4" if thrust_rk4 else "Euler"))
        plt.grid(), plt.legend(), plt.tight_layout()
        plt.savefig(sys.path[0] + "/plots/setup/integrator/benchmark/error_dt_%.4e.pdf" % dt)
        plt.close()

fig, ax = plt.subplots(figsize=(10, 6))
if only_thrust:
    con.close()
    ax.plot(fevals, mass_errors, marker="o", linewidth=1.5, label="Mass [kg]")
    ax.plot(fevals, magnitude_errors, marker="o", linewidth=1.5, label="Magnitude [N]")
else:
    ax.plot(fevals, position_errors, marker="o", linewidth=1.5, label="Position [m]")
    ax.plot(fevals, velocity_errors, marker="o", linewidth=1.5, label="Velocity [m/s]")
    if powered:
        ax.plot(fevals, mass_errors, marker="o", linewidth=1.5, label="Mass [kg]")
ax.set_xlabel("Function evaluations [-]"), ax.set_ylabel("Final error")
ax.set_xscale("log"), ax.set_yscale("log")
plt.suptitle("Using %s integrator" % ("RK 4" if thrust_rk4 else "Euler"))
plt.grid(), plt.legend(), plt.tight_layout()
plt.savefig(sys.path[0] + "/plots/setup/integrator/benchmark/dt_vs_error.pdf")