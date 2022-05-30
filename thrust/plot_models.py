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
from thrust.solid_thrust import SRM_thrust

# Specify which geometry to test
geometries_to_plot = {
    "tubular": 1,
    "rod_and_tube": 1,
    "multi_fin": 1,
    "anchor": 1,
    "spherical": 1,
    "STAR_20": 1,
    "stage_test": 1
}

print_and_show_analysis = False

def plot_geometry_and_thrust(save_path, b_s, burn_times, magnitudes, SRM_geometry, n_frames=50):
    # Resample thrust properties vs time list
    new_b_s = np.linspace(burn_times[0], burn_times[-1], n_frames)
    b_s_interp =  interpolate.interp1d(burn_times, b_s)
    mag_interp = interpolate.interp1d(burn_times, magnitudes)
    b_s = [b_s_interp(t) for t in new_b_s]
    magnitudes = [mag_interp(t) for t in new_b_s]
    burn_times = new_b_s
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
        # Add a global figure title
        fig.suptitle(str(SRM_geometry))
        fig.tight_layout()
        # Save (or show) the figure
        if save_path == "SHOW":
            plt.show()
        else:
            plt.savefig(save_path.replace(".", "-") + "%.4f.png" % burn_times[i])
            plt.close()

def to_gif(name):
    f_path = sys.path[0] + "/thrust/burn_visu/" + name
    # Convert png to gif
    os.system("convert `ls -v %s_*.png` %s.gif" % (f_path, f_path))
    # Convert gif to mov
    os.system("convert -delay 10 %s.gif %s.mov" % (f_path, f_path))
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
from models.spherical import spherical_SRM
spherical_test = spherical_SRM(R_o=0.175, R_i=0.095, run_checks=False)

STAR_20 = multi_fin_SRM(R_o=0.24, R_i=0.125, N_f=10, w_f=0.025, L_f=0.05, L=1.0, run_checks=False)

stage_test = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05, run_checks=False)
# stage_test = rod_and_tube_SRM(R_o=0.24, R_mid=0.19, R_i=0.075, L=1.05, run_checks=False)
# stage_test = anchor_SRM(R_o=0.24, R_i=0.075, N_a=3, w=0.05, r_f=0.005, delta_s=0.035, L=1.05)

geometries = [tubular_test, rod_and_tube_test, multi_fin_test, anchor_test, spherical_test, STAR_20, stage_test]

# Loop trough the geometries
for i, (name, to_plot) in enumerate(geometries_to_plot.items()):
    if to_plot:
        SRM_geometry = geometries[i]

        # Compute the thrust over time
        print("Computing thrust over time for the %s SRM geometry..." % name.replace("_", " "))
        SRM_thrust_model = SRM_thrust(SRM_geometry)
        if print_and_show_analysis:
            print("%.2f kg of propellant"%SRM_thrust_model.M_p, "%.2f kg of innert mass"%SRM_thrust_model.M_innert)
        burn_times, magnitudes, b_s, p_c_s, M_p_s = SRM_thrust_model.simulate_full_burn(compute_dep_vars=True)

        if print_and_show_analysis:
            print("Thrust of %.2f/9.75 kN for %.1f/55s" % (max(magnitudes)/1e3, burn_times[-1]))
            plt.plot(burn_times, M_p_s)
            plt.xlabel("Time [s]"), plt.ylabel("Propellant mass [kg]")
            plt.grid(), plt.tight_layout()
            plt.show()
            plot_geometry_and_thrust("SHOW", b_s, burn_times, magnitudes, SRM_geometry)

        # Plot the geometry burn over time
        print("Plotting thrust and geometry over time for the %s SRM geometry..." % name.replace("_", " "))
        plot_geometry_and_thrust(sys.path[0] + "/thrust/burn_visu/%s_" % name, b_s, burn_times, magnitudes, SRM_geometry)

        # Convert the plots to a gif
        print("Converting plots to GIF for the %s SRM geometry..." % name.replace("_", " "))
        to_gif(name)