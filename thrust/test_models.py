import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Standard imports
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

# Import solid thrust model
from thrust import solid_thrust as ST

# TODO Add main title (SRM geometry parameters)
# TODO Change burning color as a function of p_c

# Specify which geometry to test
geometries_to_plot = {
    "tubular": 1,
    "rod_and_tube": 1,
    "multi_fin": 1,
    "anchor": 1
}

def compute_thrust(SRM_geometry):
    SRM_thrust = ST.thrust(SRM_geometry)
    times = np.arange(0, 60, 0.05)
    burn_times, magnitudes, b_s, p_c_s = [], [], [], []
    for time in times:
        F_T = SRM_thrust.magnitude(time)
        magnitudes.append(F_T)
        b_s.append(SRM_thrust.b)
        burn_times.append(time)
        p_c_s.append(SRM_thrust.p_c)
        # Stop the thrust if the magnitude is 0 for the last 10 steps
        if np.sum(magnitudes[-10:]) == 0:
            break

    # Resample thrust properties vs time list
    new_burn_times = np.linspace(burn_times[0], burn_times[-1], 50)
    magnitudes = [interpolate.interp1d(burn_times, magnitudes)(t) for t in new_burn_times]
    b_s = [interpolate.interp1d(burn_times, b_s)(t) for t in new_burn_times]
    p_c_s = [interpolate.interp1d(burn_times, p_c_s)(t) for t in new_burn_times]
    burn_times = new_burn_times
    return burn_times, magnitudes, b_s, p_c_s

def plot_geometry_and_thrust(save_path, b_s, burn_times, magnitudes, SRM_geometry):
    for i, b in enumerate(b_s):
        # Create figure
        fig = plt.figure(figsize=(12,6))
        ax1 = fig.add_subplot(121, projection="polar")
        ax2 = fig.add_subplot(122)
        # Plot geometry
        SRM_geometry.plot_geometry(b, ax_in=ax1, final_b=b_s[-1])
        # Plot thrust
        ax2.plot(burn_times, np.array(magnitudes)/1e3, linewidth=1.2, color="black")
        ax2.scatter(burn_times[i], np.array(magnitudes[i])/1e3, color="black")
        ax2.plot(burn_times[:i+1], np.array(magnitudes[:i+1])/1e3, linewidth=2.0, color="black")
        ax2.set_xlabel("Time [s]"), ax2.set_ylabel("Thrust [kN]")
        ax2.grid()
        ax2.set_title("$t = %.3f$ [s]" % burn_times[i])
        fig.tight_layout()
        # Save the figure
        plt.savefig(save_path.replace(".", "-") + "%.4f.png" % burn_times[i])
        plt.close()

def to_gif(name):
    f_path = sys.path[0] + "/thrust/burn_visu/" + name
    # Duplicate first and last frames

    # Convert png to gif
    os.system("convert `ls -v %s_*.png` %s.gif" % (f_path, f_path))
    # Remove png
    os.system("rm -rf %s_*.png" % f_path)

# Define all the geometries
from models.tubular import tubular_SRM
tubular_test = tubular_SRM(R_o=1, R_i=0.25, L=3, run_checks=False)
from models.rod_and_tube import rod_and_tube_SRM
rod_and_tube_test = rod_and_tube_SRM(R_o=1, R_mid=0.65, R_i=0.25, L=3, run_checks=False)
from models.multi_fin import multi_fin_SRM
multi_fin_test = multi_fin_SRM(R_o=1, R_i=0.65, N_f=8, w_f=0.15, L_f=0.3, L=3, run_checks=False)
from models.anchor import anchor_SRM
anchor_test = anchor_SRM(R_o=1, R_i=0.25, N_a=3, w=0.2, r_f=0.05, delta_s=0.15, L=3, run_checks=False)
geometries = [tubular_test, rod_and_tube_test, multi_fin_test, anchor_test]

# Loop trough the geometries
for i, (name, to_plot) in enumerate(geometries_to_plot.items()):
    if to_plot:
        SRM_geometry = geometries[i]

        # Compute the thrust over time
        print("Computing thrust over time for the %s SRM geometry..." % name.replace("_", " "))
        burn_times, magnitudes, b_s, p_c_s = compute_thrust(SRM_geometry)

        # Plot the geometry burn over time
        print("Plotting thrust and geometry over time for the %s SRM geometry..." % name.replace("_", " "))
        plot_geometry_and_thrust(sys.path[0] + "/thrust/burn_visu/%s_" % name, b_s, burn_times, magnitudes, SRM_geometry)

        # Convert the plots to a gif
        print("Converting plots to GIF for the %s SRM geometry..." % name.replace("_", " "))
        to_gif(name)