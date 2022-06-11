# TODO:
# Plot score vs each (set of) design variables
# Plot Pareto fronts (vs relevant design variables)

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
import matplotlib.pyplot as plt
import sqlite3
import pandas as pd
import seaborn as sns

# Parameters
plot_trajectories = False
analyse_correlation = False
get_initial_population = False


# Connect to database
con = sqlite3.connect(sys.path[0]+"/optimisation/space_exploration.db")
cur = con.cursor()

if plot_trajectories:
    for dv_used in ["all", "init_angle_only", "TVC_only", "SRM_only", "*"]:
        if dv_used == "*":
            req = "SELECT id FROM solutions_multi_fin"
            line_thickness = 0.05
            alpha = 0.5
        else:
            req = "SELECT id FROM solutions_multi_fin WHERE h_a > 200e3 AND h_p > 200e3 AND h_a < 500e3 AND h_p < 500e3 AND dv_used = '%s'"%dv_used
            line_thickness = 0.5
            alpha = 1.0
        cur.execute(req)
        ids = [i[0] for i in cur.fetchall()]
        # Plot trajectories
        plt.figure(figsize=(9,5))
        for i, id in enumerate(ids):
            if i % 25 == 0:
                print("Plotting trajectory %i of %i for %s"%(i+1, len(ids), dv_used), end="\r")
            # Load simulation results
            sim_results = np.load(sys.path[0]+"/optimisation/sim_results/%i.npz"%id)
            times = sim_results["times"]
            altitudes = sim_results["dep_vars"][:,0]
            plt.plot(times[:1000]/60, altitudes[:1000]/1e3, linewidth=line_thickness, color="C%i"%(i+1), alpha=alpha)
            plt.plot(times[1001:]/60, altitudes[1001:]/1e3, linewidth=line_thickness, color="C%i"%(i+1), alpha=alpha)
            plt.plot(times[998:1002]/60, altitudes[998:1002]/1e3, linewidth=line_thickness, linestyle="dashed", color="C%i"%(i+1), alpha=alpha)
        print("Plotted %i trajectories for %s                     "%(len(ids), dv_used))
        # Prettify plot
        plt.xlabel("Time [min]")
        plt.ylabel("Altitude [km]")
        if dv_used == "*":
            plt.hlines(1e3, -1e3, 1e3, color="black", linestyle="dotted", linewidth=1.0)
            plt.ylim((-10,5500))
            plt.xlim((-2,162))
            plt.yscale("symlog", linthresh=1000, linscale=1.5)
            plt.yticks([200, 400, 600, 800, 1e3, 3e3, 5e3], ["200", "400", "600", "800", "1 $\cdot 10^3$", "3 $\cdot 10^3$", "5 $\cdot 10^3$"])
        else:
            plt.ylim((-10,502))
            plt.xlim((-2,162))
        plt.grid()
        plt.tight_layout()
        if dv_used == "*":
            plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/altitude.pdf")
            # plt.axis('off')
            # plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/altitude", transparent=True, dpi=350)
        else:
            plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/altitude_%s_sampled.pdf"%dv_used)
        plt.close()

if analyse_correlation:
    for dv_used in ["init_angle_only", "TVC_only", "SRM_only", "all"]:
        if dv_used == "init_angle_only":
            vars = ["h_a", "h_p", "angle_1", "angle_2"]
            n_rows = 2
            figsize = (4,5)
        elif dv_used == "TVC_only":
            vars = ["h_a", "h_p", "inclination"]
            n_rows = 3
            [vars.append("TVC_y_%i"%i) for i in range(1,6)]
            [vars.append("TVC_z_%i"%i) for i in range(1,6)]
            figsize = (9, 4.5)
        elif dv_used == "SRM_only":
            vars = ["h_a", "h_p", "mass"]
            n_rows = 3
            [vars.append("spherical_motor_%i"%i) for i in range(1,3)]
            [vars.append("multi_fin_motor_%i"%i) for i in range(1,7)]
            figsize = (8, 6)
        else:
            vars = ["h_a", "h_p", "mass", "angle_1", "angle_2"]
            [vars.append("TVC_y_%i"%i) for i in range(1,6)]
            [vars.append("TVC_z_%i"%i) for i in range(1,6)]
            [vars.append("spherical_motor_%i"%i) for i in range(1,3)]
            [vars.append("multi_fin_motor_%i"%i) for i in range(1,7)]
            figsize = (10, 5)
            n_rows = 3
        replace_keys = {
            "spherical_motor_1": "spherical_motor R_o",
            "spherical_motor_2": "spherical_motor R_i",
            "multi_fin_motor_1": "multi_fin_motor L",
            "multi_fin_motor_2": "multi_fin_motor R_o",
            "multi_fin_motor_3": "multi_fin_motor R_i",
            "multi_fin_motor_4": "multi_fin_motor L_f",
            "multi_fin_motor_5": "multi_fin_motor w_f",
            "multi_fin_motor_6": "multi_fin_motor N_f"
        }
    

        # Load entire database into dataframe
        req = "SELECT %s FROM solutions_multi_fin WHERE h_a_score < 1e2 AND h_p_score < 1e2 AND dv_used = '%s'"%(", ".join(vars), dv_used)
        df = pd.read_sql_query(req, con)
        # Replace dataframe column names
        df.columns = [replace_keys.get(c, c) for c in df.columns]
        # Make heatmap
        print("Making heatmap for %s..."%dv_used)
        correlations = df.corr()
        # Remove rows and columns that we don't want to compare
        correlations = correlations[:n_rows]
        correlations = correlations.iloc[:,n_rows:]
        plt.figure(figsize=figsize)
        # Plot heatmap with horizontal colorbar and xticks on top
        ax = sns.heatmap(correlations, annot=True, fmt=".2f", cmap="vlag", square=True, vmin=-1, vmax=1, cbar_kws={"orientation": "horizontal"})
        ax.xaxis.tick_top()
        ax.set_xticklabels(ax.get_xticklabels(), rotation=77.5)
        plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/heatmap_%s.pdf"%dv_used)
        plt.close()

        objectives = ["h_p_score+h_a_score"]
        vars = []
        treshold = 20
        if dv_used in ["init_angle_only", "TVC_only"]:
            operator = np.rad2deg
            if dv_used == "TVC_only":
                [vars.append("TVC_z_%i"%i) for i in range(1,6)]
                units = ["deg"]*5
            else:
                vars = ["angle_1", "angle_2"]
                units = ["deg"]*2
        elif dv_used == "SRM_only":
            treshold = 20
            objectives = ["h_p_score+h_a_score", "mass_score"]
            operator = lambda x:x
            units = ["-", "-", "m", "m", "-", "-", "-", "-"]
            [vars.append("spherical_motor_%i"%i) for i in range(1,3)]
            [vars.append("multi_fin_motor_%i"%i) for i in range(1,4)]
        else:
            continue

        # Load entire database into dataframe
        req = "SELECT %s, %s FROM solutions_multi_fin WHERE h_a_score < 1e2 AND h_p_score < 1e2 AND dv_used = '%s'"%(", ".join(objectives), ", ".join(vars), dv_used)
        # if dv_used == "SRM_only":
        #     req = "SELECT %s, %s FROM solutions_multi_fin WHERE h_a_score < 1e2 AND h_p_score < 1e2 AND dv_used = '%s'"%(", ".join(objectives), ", ".join(vars), dv_used)
        df = pd.read_sql_query(req, con)
        # Extract objective columns
        df_objectives = df[objectives]
        # Replace dataframe column names
        df.columns = [replace_keys.get(c, c) for c in df.columns]
        # Replace name in vars
        vars = [replace_keys.get(c, c) for c in vars]
        # Plot objective score vs design variables
        print("Plotting objective score vs design variables for %s..."%dv_used)
        for objective in objectives:
            plt.figure(figsize=(9, 5))
            for i, var in enumerate(vars):
                f = 10 if "N_f" in var else 1
                # Index of objectives below treshold
                idx = df_objectives[objective] <= treshold
                # Variables of objectives below treshold
                x = operator(df[var][idx])
                # Maximum and minimum variable value
                x_max, x_min = np.max(x), np.min(x)
                print("For %s below treshold, %s has range %.2f to %.2f" % (objective, var, x_min, x_max))
                plt.plot(operator(df[var])/f, df_objectives[objective], "o", label="%s [%s]"%("%s/%i"%(var, f) if f != 1 else var, units[i]), markersize=1.5, alpha=0.5)
            plt.grid()
            plt.legend()
            plt.yscale("log")
            plt.xlabel("Design variable value")
            ylims = plt.ylim()
            plt.ylim((1e-1, ylims[1]))
            plt.ylabel(objective)
            plt.tight_layout()
            plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/%s_vs_%s_dvs.pdf"%(objective, dv_used))

if get_initial_population:
    N = 48
    samples_weights = {
        "init_angle_only": 2,
        "TVC_only": 1,
        "SRM_only": 3,
        "all": 1
    }
    sum_weights = sum(samples_weights.values())
    samples_weights = {k: v/sum_weights for k, v in samples_weights.items()}
    plt.figure(figsize=(9,5))
    i_tot = 0
    for dv_used, fraction in samples_weights.items():
        sim_num = int(N*fraction)
        if dv_used == "all":
            sim_num = N - i_tot
        req = "SELECT id FROM solutions_multi_fin WHERE dv_used = '%s' ORDER BY h_a_score+h_p_score+mass_score ASC LIMIT %i"%(dv_used, sim_num)
        cur.execute(req)
        ids = [i[0] for i in cur.fetchall()]
        for i, id in enumerate(ids):
            # Load simulation results
            sim_results = np.load(sys.path[0]+"/optimisation/sim_results/%i.npz"%id)
            times = sim_results["times"]
            altitudes = sim_results["dep_vars"][:,0]
            plt.plot(times[:1000]/60, altitudes[:1000]/1e3, linewidth=0.15, color="C%i"%(i+1))
            plt.plot(times[1001:]/60, altitudes[1001:]/1e3, linewidth=0.15, color="C%i"%(i+1))
            plt.plot(times[998:1002]/60, altitudes[998:1002]/1e3, linewidth=0.15, linestyle="dashed", color="C%i"%(i+1))
        print("Plotted %i trajectories for initial population from %s"%(len(ids), dv_used))
        i_tot += len(ids)
    print("Plotted %i trajectories for initial population"%i_tot)
    # Prettify plot
    plt.xlabel("Time [min]")
    plt.ylabel("Altitude [km]")
    plt.grid()
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/initial_population.pdf")

# Close the database connection
con.close()