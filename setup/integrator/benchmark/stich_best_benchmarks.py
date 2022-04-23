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
import glob
import numpy as np
from matplotlib import pyplot as plt


times_stiched, states_stiched, dep_vars_stiched = [], [], []
f_evals_stiched = 0

# Loop trough all stage and phases of the ascent
for stage in [1, 2]:
    for powered in [True, False]:
        # Load the corresponding benchmark result file
        f_name = "/data/best_integrator_dt/%i_%s_dt_*.npz" % (stage, "V" if powered else "X")
        try:
            last_state_fname = glob.glob(sys.path[0]+f_name)[0]
        except IndexError:
            print("No files found for stage %i, %s" % (stage, "powered" if powered else "unpowered"))
            continue
        saved_results = np.load(last_state_fname)

        # Extract the times and states
        times = saved_results["times"]
        states = saved_results["states"]
        dep_vars = saved_results["dep_vars"]
        f_evals = saved_results["f_evals"]
        altitudes = dep_vars[:, 0]

        # Save the final mass if thrust was used
        if powered:
            final_mass = states[-1,-1]
        # Add the constant mass to the states if no thrust was used
        else:
            states = np.concatenate([states, np.ones((states.shape[0],1),dtype=states.dtype)*final_mass], axis=1)

        # Print the initial mass of the current ascent part
        print("Initial mass for stage %i %s [kg] =" % (stage, "powered" if powered else "unpowered"), states[0,-1])

        # Save the times and states
        times_stiched.append(times)
        dep_vars_stiched.append(altitudes)
        states_stiched.append(states)
        f_evals_stiched += f_evals


# Reshape the arrays
times_stiched = np.asarray(times_stiched).flatten()
dep_vars_stiched = np.asarray(dep_vars_stiched).reshape((len(times_stiched), 1))
states_stiched = np.asarray(states_stiched).reshape((len(times_stiched), -1))

# Save the full benchmark
np.savez(sys.path[0]+"/data/best_integrator_dt/full_benchmark.npz", times=times_stiched, states=states_stiched, dep_vars=dep_vars_stiched, f_evals=f_evals_stiched)

# Plot benchmark results
fig = plt.figure(figsize=(9, 5))
gs = fig.add_gridspec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
# Only plot mass up to 13min
idx_crop = np.where(times_stiched >= 13*60)[0][0]
ax1.plot(times_stiched/60, dep_vars_stiched[:,0]/1e3, label="Benchmark")
ax2.plot(times_stiched[:idx_crop]/60, states_stiched[:idx_crop,-1], label="Benchmark")

# Load and plot ascent test results
f_name = "/data/MAV_ascent_test.npz"
saved_results = np.load(sys.path[0]+f_name)
times = saved_results["times"]
states = saved_results["states"]
dep_vars = saved_results["dep_vars"]
altitudes = dep_vars[:, 1]
masses = states[:, -1]
idx_crop = np.where(times >= 13*60)[0][0]
ax1.plot(times/60, altitudes/1e3, label="Variable step test", linestyle="dotted")
ax2.plot(times[:idx_crop]/60, masses[:idx_crop], label="Variable step test", linestyle="dotted")

# Prettify plot
ax1.set_xlabel("Time [min]")
ax1.set_ylabel("Altitude [km]")
ax2.set_xlabel("Time [min]")
ax2.set_ylabel("Mass [kg]")
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.show()