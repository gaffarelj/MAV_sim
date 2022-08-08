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
import sqlite3
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Custom imports
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.tubular import tubular_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.anchor import anchor_SRM
from optimisation.ascent_problem import MAV_problem

# Tudatpy imports
from tudatpy import util

N = 200
run_high_accuracy_sims = False
compare_accuracy = False

if __name__ == "__main__":
    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
    cur = con.cursor()
    fig1 = plt.figure(figsize=(9,3.5))
    fig2 = plt.figure(figsize=(9,3.5))
    i_plt = 0
    df_diffs = pd.DataFrame(columns=["normal_id", "accurate_id", "diff_pos", "diff_vel"])

    for SRM_type in ["multi_fin", "anchor", "tubular", "rod_and_tube"]:

        # Get N results with h_a_score+h_p_score = 0 and mass_score < 1.0
        req = "SELECT * FROM solutions_%s WHERE angle_1 IS NOT NULL and h_a_score IS NOT NULL AND h_a_score + h_p_score = 0 AND mass_score < 1.0 ORDER BY mass_score ASC LIMIT %i"%(SRM_type, N//4)
        df = pd.read_sql_query(req, con)

        if run_high_accuracy_sims:
            dv_s = []
            n_dv = None
            for i_row in range(len(df)):
                if SRM_type == "multi_fin":
                    dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["multi_fin_motor_1"].iloc[i_row], df["multi_fin_motor_2"].iloc[i_row], df["multi_fin_motor_3"].iloc[i_row], df["multi_fin_motor_4"].iloc[i_row], df["multi_fin_motor_5"].iloc[i_row], df["multi_fin_motor_6"].iloc[i_row]]
                elif SRM_type == "tubular":
                    dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["tubular_motor_1"].iloc[i_row], df["tubular_motor_2"].iloc[i_row], df["tubular_motor_3"].iloc[i_row]]
                elif SRM_type == "rod_and_tube":
                    dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["rod_and_tube_motor_1"].iloc[i_row], df["rod_and_tube_motor_2"].iloc[i_row], df["rod_and_tube_motor_3"].iloc[i_row], df["rod_and_tube_motor_4"].iloc[i_row]]
                elif SRM_type == "anchor":
                    dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["anchor_motor_1"].iloc[i_row], df["anchor_motor_2"].iloc[i_row], df["anchor_motor_3"].iloc[i_row], df["anchor_motor_4"].iloc[i_row], df["anchor_motor_5"].iloc[i_row], df["anchor_motor_6"].iloc[i_row], df["anchor_motor_7"].iloc[i_row]]
                dv_s.extend(dv)
                if n_dv is None:
                    n_dv = len(dv)

            if SRM_type == "multi_fin":
                SRM_class = multi_fin_SRM
            elif SRM_type == "anchor":
                SRM_class = anchor_SRM
            elif SRM_type == "tubular":
                SRM_class = tubular_SRM
            elif SRM_type == "rod_and_tube":
                SRM_class = rod_and_tube_SRM

            problem = MAV_problem([[1]*n_dv]*2, SRM_class)
            # problem.batch_fitness(dv_s, plot_geo=True, plot_thrust=True, better_accuracy=True, save_to_db="extra_accurate")
            problem.batch_fitness(dv_s, plot_geo=False, plot_thrust=False, better_accuracy=True, save_to_db="extra_accurate")

        elif compare_accuracy:
            for i_row in range(len(df)):
                req = "SELECT * FROM solutions_%s WHERE mass_score = %s AND dv_used = 'extra_accurate'"%(SRM_type, df["mass_score"].iloc[i_row])
                df_accurate = pd.read_sql_query(req, con)
                # Load results and plot difference
                try:
                    results_normal = np.load(sys.path[0]+"/optimisation/sim_results/%s/%i.npz"%(SRM_type, df["id"].iloc[i_row]))
                    results_accurate = np.load(sys.path[0]+"/optimisation/sim_results/%s/%i.npz"%(SRM_type, df_accurate["id"].iloc[0]))
                except IndexError:
                    continue
                dict_normal = {t: state for t, state in zip(results_normal["times"], results_normal["states"])}
                dict_accurate = {t: state for t, state in zip(results_accurate["times"], results_accurate["states"])}
                dict_diff = util.compare_results(dict_normal, dict_accurate, np.linspace(1, 159, 500))
                diff_states = util.result2array(dict_diff)
                diff_pos, diff_vel = np.linalg.norm(diff_states[:,1:4], axis=1), np.linalg.norm(diff_states[:,4:7], axis=1)
                plt.figure(fig1)
                plt.plot(diff_states[:,0], diff_pos, linewidth=0.5)
                plt.figure(fig2)
                plt.plot(diff_states[:,0], diff_vel, linewidth=0.5)
                df_diffs.loc[i_plt] = [
                    df["id"].iloc[i_row],
                    df_accurate["id"].iloc[0],
                    diff_pos[-1],
                    diff_vel[-1]
                ]
                i_plt += 1

    if compare_accuracy:
        plt.figure(fig1)
        plt.ylabel("Absolute difference in position [m]"), plt.xlabel("Time [min]")
        plt.yscale("log"), plt.grid(), plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimums_accuracy_comparison_pos.pdf")
        plt.figure(fig2)
        plt.ylabel("Absolute difference in velocity [m/s]"), plt.xlabel("Time [min]")
        plt.yscale("log"), plt.grid(), plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimums_accuracy_comparison_vel.pdf")
        print(df_diffs.describe())
        print(i_plt)

    con.close()
# TODO: Ideally, plot the q0.25, q0.5, q0.75, (compute percentiles and take DV that was actually run) etc that represent the most common geometries for both stage (plot geo in a big subplot)