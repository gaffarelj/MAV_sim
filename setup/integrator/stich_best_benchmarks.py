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

# Load and plot benchmark results from each ascent part
fig = plt.figure(figsize=(14, 5))
gs = fig.add_gridspec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
stage_1_mass_hist, stage_2_mass_hist = [], []
times_stiched = []
states_stiched = []
dep_vars_stiched = []
for stage in [1, 2]:
    for powered in [True, False]:
        f_name = "/data/best_integrator_dt/%i_%s_dt_*.npz" % (stage, "V" if powered else "X")
        try:
            last_state_fname = glob.glob(sys.path[0]+f_name)[0]
        except IndexError:
            print("No files found for stage %i, %s" % (stage, "powered" if powered else "unpowered"))
            continue
        saved_results = np.load(last_state_fname)
        times = saved_results["times"]
        states = saved_results["states"]
        dep_vars = saved_results["dep_vars"]
        altitudes = dep_vars[:, 0]

        # Plot altitude as a function of time
        ax1.plot(times, altitudes/1e3, label="Stage %i, %s" % (stage, "powered" if powered else "unpowered"))
        times_stiched.append(times)
        dep_vars_stiched.append(altitudes)
        if stage == 2 and powered:
            m_times = [stage_1_mass_hist[0], times[0]]
            masses = [stage_1_mass_hist[1], stage_1_mass_hist[1]]
            ax2.plot(m_times, masses, label="Stage 1, unpowered")
        if powered:
            masses = states[:, -1]
            if stage == 1:
                ax2.plot(times, masses, label="Stage 1, powered")
            else:
                ax2.plot([times[0], *times], [stage_1_mass_hist[1], *masses], label="Stage 2, powered")
        elif stage == 2:
            ax2.plot([times[0], times[-1]], [stage_2_mass_hist[0], stage_2_mass_hist[0]], label="Stage 2, unpowered")
        if stage == 1 and not powered:
            stage_1_mass_hist = [times[0], masses[-1]]
        elif stage == 2 and powered:
            stage_2_mass_hist = [masses[-1]]

times_stiched = np.asarray(times_stiched).flatten()
dep_vars_stiched = np.asarray(dep_vars_stiched).reshape((len(times_stiched), 1))

# print(times_stiched, dep_vars_stiched)

# input()
# ax1.plot(times_stiched, dep_vars_stiched[:,0]/1e3, label="Stiched")

# Load and plot ascent test results
f_name = "/data/MAV_ascent_test.npz"
saved_results = np.load(sys.path[0]+f_name)
times = saved_results["times"]
states = saved_results["states"]
dep_vars = saved_results["dep_vars"]
altitudes = dep_vars[:, 1]
masses = states[:, -1]
ax1.plot(times, altitudes/1e3, label="Ascent test", linestyle="dotted")
ax2.plot(times, masses, label="Ascent test", linestyle="dotted")

ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Altitude [km]")
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("Mass [kg]")
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.show()