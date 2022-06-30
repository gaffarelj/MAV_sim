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
import re

# Tudatpy imports
from tudatpy import plotting

# Custom imports
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM

# Connect to the database
con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
cur = con.cursor()

def get_df_gens(df):
    # Extract generation number from dv_used column, and replace column name
    df["dv_used"] = df["dv_used"].str.split("_").str[-1].astype(int)
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "multi_fin_motor_1": "multi_fin_motor_L", "multi_fin_motor_2": "multi_fin_motor_R_o", "multi_fin_motor_3": "multi_fin_motor_R_i", "multi_fin_motor_4": "multi_fin_motor_L_f", "multi_fin_motor_5": "multi_fin_motor_w_f", "multi_fin_motor_6": "multi_fin_motor_N_f"}
    df.rename(columns=replace_cols, inplace=True)
    return df



for seed in [42, 13, 123, 846, 579]:
    # Load optimisation results
    req = "SELECT * FROM solutions_multi_fin WHERE dv_used LIKE 'opti_%s_%%' AND h_p_score IS NOT NULL AND h_p_score < 5e2 AND h_a_score < 5e2" % seed
    df = pd.read_sql_query(req, con)
    df = get_df_gens(df)

    des_vars_legends = [("$\\theta_1$", "$\\theta_2$"), ("$TVC_{z,1}$", "$TVC_{z,2}$", "$TVC_{z,3}$", "$TVC_{z,4}$", "$TVC_{z,5}$"), ("$R_o$", "$R_i$"), ("$L$", "$R_o$", "$R_i$", "$L_f$", "$w_f$", "$N_f$")]
    des_vars_legends_flat = [item for tup in des_vars_legends for item in tup]
    des_vars_names = ["angle_1", "angle_2", "TVC_z_1", "TVC_z_2", "TVC_z_3", "TVC_z_4", "TVC_z_5", "spherical_motor_R_o", "spherical_motor_R_i", "multi_fin_motor_L", "multi_fin_motor_R_o", "multi_fin_motor_R_i", "multi_fin_motor_L_f", "multi_fin_motor_w_f", "multi_fin_motor_N_f"]

    # Define the range in which the design variables can vary
    launch_angle_1_range = np.deg2rad([47.5, 60])
    launch_angle_2_range = np.deg2rad([70, 110])
    TVC_range = np.deg2rad([-5, 5])
    N_TVC_nodes = 5
    spherical_SRM_range = [[0.3, 0.2], [1.0, 0.9]]
    multi_fin_SRM_range = [[0.3, 0.1, 0.2, 0.25, 0.35, 3], [1.25, 0.285, 0.9, 0.75, 0.9, 20]]
    design_var_range = (
        [launch_angle_1_range[0], launch_angle_2_range[0], *[TVC_range[0]]*N_TVC_nodes, *spherical_SRM_range[0], *multi_fin_SRM_range[0]],
        [launch_angle_1_range[1], launch_angle_2_range[1], *[TVC_range[1]]*N_TVC_nodes, *spherical_SRM_range[1], *multi_fin_SRM_range[1]]
    )

    def rescale_design_vars(df):
        # Rescale values to the range of the design variables
        for i, tag in enumerate(des_vars_names):
            df[tag] = (df[tag] - design_var_range[0][i])/(design_var_range[1][i] - design_var_range[0][i])
        return df

    # Plot the history of the objective scores
    mean_scores = df.groupby("generation").mean()
    min_scores = df.groupby("generation").min()
    plt.figure(figsize=(9, 5))
    plt.plot(mean_scores["h_p_score"], label="$h_p$ score (mean)", color="C0", linestyle="solid")
    plt.plot(min_scores["h_p_score"], label="$h_p$ score (min)", color="C0", linestyle="dashed")
    plt.plot(mean_scores["h_a_score"], label="$h_a$ score (mean)", color="C1", linestyle="solid")
    plt.plot(min_scores["h_a_score"], label="$h_a$ score (min)", color="C1", linestyle="dotted")
    plt.plot(mean_scores["mass_score"], label="mass score (mean)", color="C2", linestyle="solid")
    plt.plot(min_scores["mass_score"], label="mass score (min)", color="C2", linestyle="dashed")
    plt.grid()
    plt.legend(loc="upper left")
    plt.xlabel("Generation [-]")
    plt.ylabel("Objective score [-]")
    plt.ylim([-0.5, 10])
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/results/history_mm_%i.pdf" % seed)

    # Index of points with h_p_score + h_a_score < 1.5
    df["h_comb_scores"] = df["h_p_score"] + df["h_a_score"]
    idx_h_score = df["h_comb_scores"].values < 1.5
    # Index of points with mass_score < 1.0
    idx_mass_score = df["mass_score"].values < 1.0
    combined_idxs = idx_h_score & idx_mass_score

    # Plot the Pareto fronts
    fig, ax = plotting.pareto_front(
        x_objective=df["h_a_score"].loc[combined_idxs]+df["h_p_score"].loc[combined_idxs],
        y_objective=df["mass_score"].loc[combined_idxs],
        x_label="Altitude score (apoapsis + periapsis)",
        y_label="Mass score",
        alpha=0.65
    )
    # Save the plot
    plt.savefig(sys.path[0]+"/plots/optimisation/results/Pareto_%i.pdf" % seed)
    plt.close()

    df = rescale_design_vars(df)

    # Plot history of design variables
    std_des_vars = df.groupby("generation").std()
    plt.figure(figsize=(9, 5))
    for i_dv, des_var in enumerate(des_vars_names):
        try:
            ls = "dashed" if i_dv < 10 else "dotted"
            plt.plot(std_des_vars[des_var], label=des_vars_legends_flat[i_dv], color="C%d" % i_dv, linestyle=ls)
        except KeyError:
            pass
    plt.grid()
    plt.legend(loc="upper center", ncol=5)
    plt.xlabel("Generation [-]")
    plt.ylabel("$\sigma$ of design variable scaled in their range [-]")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/results/history_std_%i.pdf" % seed)

## Get optimum values
# Select where h_comb_score = 0 and mass_score is min
req = "SELECT * FROM solutions_multi_fin WHERE dv_used LIKE 'opti_%%' AND h_p_score+h_a_score = 0 ORDER BY mass_score ASC LIMIT 3"
df_opt = pd.read_sql(req, con=con)
# Order df_opt by mass_score
df_opt = df_opt.sort_values(by="mass_score")
print(df_opt.columns)
for i_row in range(df_opt.shape[0]):
    # Get the design variables
    print("Initial angles [deg]: %.2f, %.2f" % tuple([np.rad2deg(df_opt.iloc[i_row][col_n]) for col_n in ["angle_1", "angle_2"]]))
    print("TVC angles: %.2f, %.2f, %.2f, %.2f, %.2f" % tuple([np.rad2deg(df_opt.iloc[i_row]["TVC_z_%i"%i]) for i in range(1, 6)]))
    L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = [df_opt.iloc[i_row]["multi_fin_motor_%i"%i] for i in range(1, 7)]
    R_i_1 = R_i_frac_1 * R_o_1
    L_f = L_f_frac * R_i_1
    w_f = w_f_frac*2*np.pi*(R_i_1-L_f)/N_f
    # Make sure that fins will burn before outer tube
    if w_f/2 > R_o_1 - R_i_1:
        w_f = 2*(R_o_1 - R_i_1)-1e-7
    R_o_2_frac, R_i_2_frac = [df_opt.iloc[i_row]["spherical_motor_%i"%i] for i in range(1, 3)]
    R_o_2 = R_o_2_frac * R_o_1
    R_i_2 = R_i_2_frac * R_o_2
    SRM_2_model = spherical_SRM(R_o_2, R_i_2)
    SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
    SRM_str = str(SRM_1_model).replace("\n", "").replace("$", " ")
    SRM_str = re.sub(r"\s+", " ", SRM_str)
    print("First stage motor:", SRM_str)
    print("First stage motor mass [kg]: %.2f" % (SRM_1_model.get_V_p()*1854.5))
    print("First stage motor burn time [s]: %.2f" % df_opt.iloc[i_row]["t_b_1"])
    SRM_str = str(SRM_2_model).replace("\n", "").replace("$", " ")
    SRM_str = re.sub(r"\s+", " ", SRM_str)
    print("Second stage motor:", SRM_str)
    print("Second stage motor mass [kg]: %.2f" % (SRM_2_model.get_V_p()*1854.5))
    print("Second stage motor burn time [s]: %.2f" % df_opt.iloc[i_row]["t_b_2"])

    print()

con.close()


# TODO: plot group of DV for objectives = 0 for h; < 0.95 for mass, to show a good (and normally robust) design. Make one such plot per set of DV