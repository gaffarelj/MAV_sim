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

# Connect to the database
con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
cur = con.cursor()

# Parameters
analyse_history = False
analyse_algorithm = True

def plot_vert_dots(df, x_tag, y_tags, fname, legends=None, y_label="Score", x_label="Generation", figsize=(9,5), within_one=False):
    plt.figure(figsize=figsize)
    markers = ["v", "^", "P", ">", "<", "X", "D", "x", "s", "p", "*", "h", "H", "o", "d", "+", "8"]
    if legends is None:
        legends = y_tags
    for i, y_tag in enumerate(y_tags):
        indent = ((i+1)/(len(y_tags)+1)-0.5)
        if within_one:
            plt.scatter(df[x_tag]+indent, df[y_tag], s=3.0, alpha=0.45, color="C%i"%i, marker=markers[i], label=legends[i])
        else:
            plt.scatter(df[x_tag]+indent, df[y_tag], s=3.0, alpha=0.45, color="C%i"%i, marker=markers[i])
            plt.scatter(np.unique(df[x_tag])+indent, df.groupby(x_tag).min()[y_tag], s=20.0, color="C%i"%i, marker=markers[i], label=legends[i])

    # Prettify
    plt.ylabel(y_label)
    if within_one:
        plt.ylim(-0.1, 1.1)
    else:
        plt.yscale("symlog", linthresh=0.5)
    # Set vertical line to separate generations
    for x in np.arange(0.5, len(np.unique(df[x_tag]))-0.4, 1):
        plt.axvline(x, color="k", alpha=0.25, linestyle="--", linewidth=0.5)
    plt.grid(axis="y")
    plt.legend(ncol=8, loc="upper center")

    plt.xlabel(x_label)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

def get_df_gens(df):
    # Extract generation number from dv_used column, and replace column name
    df["dv_used"] = df["dv_used"].str.split("_").str[-1].astype(int)
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "multi_fin_motor_1": "multi_fin_motor_L", "multi_fin_motor_2": "multi_fin_motor_R_o", "multi_fin_motor_3": "multi_fin_motor_R_i", "multi_fin_motor_4": "multi_fin_motor_L_f", "multi_fin_motor_5": "multi_fin_motor_w_f", "multi_fin_motor_6": "multi_fin_motor_N_f"}
    df.rename(columns=replace_cols, inplace=True)
    return df

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

des_vars_names = ["angle_1", "angle_2", "TVC_z_1", "TVC_z_2", "TVC_z_3", "TVC_z_4", "TVC_z_5", "spherical_motor_R_o", "spherical_motor_R_i", "multi_fin_motor_L", "multi_fin_motor_R_o", "multi_fin_motor_R_i", "multi_fin_motor_L_f", "multi_fin_motor_w_f", "multi_fin_motor_N_f"]

def rescale_design_vars(df):
    # Rescale values to the range of the design variables
    for i, tag in enumerate(des_vars_names):
        df[tag] = (df[tag] - design_var_range[0][i])/(design_var_range[1][i] - design_var_range[0][i])
    return df

des_vars_fnames = ["angles", "TVC", "spherical_SRM", "multi_fin_SRM"]
des_vars_idxs = [(0, 2), (2, 7), (7, 9), (9, 15)]
des_vars_legends = [("$\\theta_1$", "$\\theta_2$"), ("$TVC_{z,1}$", "$TVC_{z,2}$", "$TVC_{z,3}$", "$TVC_{z,4}$", "$TVC_{z,5}$"), ("$R_o$", "$R_i$"), ("$L$", "$R_o$", "$R_i$", "$L_f$", "$w_f$", "$N_f$")]

if analyse_history:
    algorithms = ["NSGA2", "MACO", "NSPSO"]
    # Get all results for each algorithm
    for algorithm in algorithms:
        print("Analysing algorithm %s"%algorithm)

        # Load all solutions results into a dataframe
        req = "SELECT * FROM solutions_multi_fin WHERE dv_used LIKE '%s_%%_tuning_%%' AND h_p_score < 5e2 AND h_a_score < 5e2" % algorithm
        df = pd.read_sql_query(req, con)
        df = get_df_gens(df)

        # Plot the history of the objective scores
        plot_vert_dots(
            df,
            "generation",
            ["h_p_score", "h_a_score", "mass_score"],
            sys.path[0]+"/plots/optimisation/tuning/%s_history.pdf" % algorithm,
            ["$h_p$", "$h_a$", "mass"],
            "Objective score [-]"
        )

        df = rescale_design_vars(df)

        # Plot the design variable history in group
        for des_vars_fname, des_vars_idx, des_vars_legend in zip(des_vars_fnames, des_vars_idxs, des_vars_legends):
            plot_vert_dots(
                df,
                "generation",
                des_vars_names[des_vars_idx[0]:des_vars_idx[1]],
                sys.path[0]+"/plots/optimisation/tuning/%s_%s_dv.pdf" % (algorithm, des_vars_fname),
                des_vars_legend,
                "Design variable value in its range [-]",
                within_one=True
            )
    
if analyse_algorithm:    
    params_base_values = {
        "omega": 0.6,
        "c1": 0.01,
        "c2": 0.5,
        "chi": 0.5,
        "v_coeff": 0.5,
        "leader_selection_range": 2,
        "diversity_mechanism": "crowding distance"
    }
    parameter_names = list(params_base_values.keys())

    # Get all results for each parameter
    for parameter_name in parameter_names:
        print("Analysing parameter %s"%parameter_name)
        # Get parameter values
        req = "SELECT DISTINCT(dv_used) FROM solutions_multi_fin WHERE dv_used LIKE '%s_%%_pt_%%' AND h_p_score IS NOT NULL" % parameter_name
        res = cur.execute(req).fetchall()
        parameter_values = list(set([float(r[0].split("_")[1]) for r in res]))
        # Get all match for each parameter value
        for parameter_value in parameter_values:
            if parameter_value == params_base_values[parameter_name]:
                req = "SELECT * FROM solutions_multi_fin WHERE dv_used LIKE 'NSPSO_%_tuning_%'"
                print("Parameter value %s is the benchmark"%parameter_value)
            else:
                req = "SELECT * FROM solutions_multi_fin WHERE dv_used LIKE '%s_%s_%%_pt_%%' AND h_p_score IS NOT NULL AND h_p_score < 5e2 AND h_a_score < 5e2" % (parameter_name, parameter_value)
            df = pd.read_sql_query(req, con)
            df = get_df_gens(df)

            # Create param_name directory if it doesn't exist
            if not os.path.exists(sys.path[0]+"/plots/optimisation/tuning/%s" % parameter_name):
                os.makedirs(sys.path[0]+"/plots/optimisation/tuning/%s" % parameter_name)

            # Plot the history of the objective scores
            plot_vert_dots(
                df,
                "generation",
                ["h_p_score", "h_a_score", "mass_score"],
                sys.path[0]+"/plots/optimisation/tuning/%s/%s_history.pdf" % (parameter_name, parameter_value),
                ["$h_p$", "$h_a$", "mass"],
                "Objective score [-]"
            )

            df = rescale_design_vars(df)

            # Plot the design variable history in group
            for des_vars_fname, des_vars_idx, des_vars_legend in zip(des_vars_fnames, des_vars_idxs, des_vars_legends):
                plot_vert_dots(
                    df,
                    "generation",
                    des_vars_names[des_vars_idx[0]:des_vars_idx[1]],
                    sys.path[0]+"/plots/optimisation/tuning/%s/%s_%s_dv.pdf" % (parameter_name, des_vars_fname, parameter_value),
                    des_vars_legend,
                    "Design variable value in its range [-]",
                    within_one=True
                )
            
            # # Get actual DVs where None in df (was already run, sim was skipped)
            # null_idxs = list(df[df.angle_1.isnull()].index)
            # reqs = []
            # tmp_idxs = []
            # while len(null_idxs) > 0:
            #     null_idx = null_idxs.pop(0)
            #     tmp_idxs.append(null_idx)
            #     null_dvs = [df.loc[null_idx, "h_p_score"], df.loc[null_idx, "h_a_score"], df.loc[null_idx, "mass_score"]]
            #     reqs.append("SELECT * FROM solutions_multi_fin WHERE h_p_score = %s AND h_a_score = %s AND mass_score = %s AND angle_1 IS NOT NULL" % (null_dvs[0], null_dvs[1], null_dvs[2]))
            #     if len(reqs) == 495 or len(null_idxs) == 0:
            #         print("Running %i SQL requests..." % len(reqs))
            #         rep_df = pd.read_sql_query(" UNION ".join(reqs), con)
            #         for i, null_idx in enumerate(tmp_idxs):
            #             # Replace column values by request results
            #             gen_idx = df.loc[null_idx, "generation"]
            #             df.loc[null_idx] = rep_df.iloc[i]
            #             df.loc[null_idx, "generation"] = gen_idx
            #         reqs = []
            #         tmp_idxs = []

            # # if len(null_idxs) > 0:
            # #     reqs = []
            # #     for null_idx in null_idxs:
            # #         null_dvs = [df.loc[null_idx, "h_p_score"], df.loc[null_idx, "h_a_score"], df.loc[null_idx, "mass_score"]]
            # #         # Get generation number
            # #         reqs.append("SELECT * FROM solutions_multi_fin WHERE h_p_score = %s AND h_a_score = %s AND mass_score = %s AND angle_1 IS NOT NULL" % (null_dvs[0], null_dvs[1], null_dvs[2]))
            # #     rep_df = pd.read_sql_query(" UNION ".join(reqs), con)
            # #     for i, null_idx in enumerate(null_idxs):
            # #         # Replace column values by request results
            # #         gen_idx = df.loc[null_idx, "generation"]
            # #         df.loc[null_idx] = rep_df.iloc[i]
            # #         df.loc[null_idx, "generation"] = gen_idx
            
            # # print(df), input()


con.close()