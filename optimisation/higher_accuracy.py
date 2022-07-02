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
from optimisation.ascent_problem import MAV_problem

# Tudatpy imports
from tudatpy import util

N = 100
run_high_accuracy_sims = False
compare_accuracy = False

if __name__ == "__main__":
    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
    cur = con.cursor()

    # Get N results with h_a_score+h_p_score = 0 and mass_score < 1.0
    req = "SELECT * FROM solutions_multi_fin WHERE angle_1 IS NOT NULL and h_a_score IS NOT NULL AND h_a_score + h_p_score = 0 AND mass_score < 1.0 ORDER BY mass_score ASC LIMIT %i"%N
    df = pd.read_sql_query(req, con)

    if run_high_accuracy_sims:
        dv_s = []
        for i_row in range(len(df)):
            dv = default_dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["multi_fin_motor_1"].iloc[i_row], df["multi_fin_motor_2"].iloc[i_row], df["multi_fin_motor_3"].iloc[i_row], df["multi_fin_motor_4"].iloc[i_row], df["multi_fin_motor_5"].iloc[i_row], df["multi_fin_motor_6"].iloc[i_row]]
            dv_s.extend(dv)

        problem = MAV_problem([[1]*15]*2, multi_fin_SRM)
        # problem.batch_fitness(dv_s, plot_geo=True, plot_thrust=True, better_accuracy=True, save_to_db="extra_accurate")
        problem.batch_fitness(dv_s, plot_geo=False, plot_thrust=False, better_accuracy=True, save_to_db="extra_accurate")

    elif compare_accuracy:
        df_diffs = pd.DataFrame(columns=["normal_id", "accurate_id", "diff_pos", "diff_vel"])
        fig1 = plt.figure(figsize=(9,5))
        fig2 = plt.figure(figsize=(9,5))
        for i_row in range(len(df)):
            req = "SELECT * FROM solutions_multi_fin WHERE mass_score = %s AND dv_used = 'extra_accurate'"%(df["mass_score"].iloc[i_row])
            df_accurate = pd.read_sql_query(req, con)
            # Load results and plot difference
            results_normal = np.load(sys.path[0]+"/optimisation/sim_results/%i.npz"%(df["id"].iloc[i_row]))
            results_accurate = np.load(sys.path[0]+"/optimisation/sim_results/%i.npz"%(df_accurate["id"].iloc[0]))
            dict_normal = {t: state for t, state in zip(results_normal["times"], results_normal["states"])}
            dict_accurate = {t: state for t, state in zip(results_accurate["times"], results_accurate["states"])}
            dict_diff = util.compare_results(dict_normal, dict_accurate, np.linspace(1, 159, 500))
            diff_states = util.result2array(dict_diff)
            diff_pos, diff_vel = np.linalg.norm(diff_states[:,1:4], axis=1), np.linalg.norm(diff_states[:,4:7], axis=1)
            plt.figure(fig1)
            plt.plot(diff_states[:,0], diff_pos, linewidth=0.5)
            plt.figure(fig2)
            plt.plot(diff_states[:,0], diff_vel, linewidth=0.5)
            df_diffs.loc[i_row] = [
                df["id"].iloc[i_row],
                df_accurate["id"].iloc[0],
                diff_pos[-1],
                diff_vel[-1]
            ]
        plt.figure(fig1)
        plt.ylabel("Absolute difference in position [m]"), plt.xlabel("Time [min]")
        plt.yscale("log"), plt.grid(), plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimums_accuracy_comparison_pos.pdf")
        plt.figure(fig2)
        plt.ylabel("Absolute difference in velocity [m/s]"), plt.xlabel("Time [min]")
        plt.yscale("log"), plt.grid(), plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimums_accuracy_comparison_vel.pdf")
        print(df_diffs.describe())

    con.close()
# TODO: Ideally, plot the q0.25, q0.5, q0.75, (compute percentiles and take DV that was actually run) etc that represent the most common geometries for both stage (plot geo in a big subplot)