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
plot_trajectories = True
analyse_correlation = False
get_initial_population = False
pareto_fronts = False

# Tudatpy imports
from tudatpy import plotting


# Connect to database
con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db")
cur = con.cursor()

if plot_trajectories:
    nice_png_fig = True
    all_dvs_types = ["*"] if nice_png_fig else ["all", "init_angle_only", "TVC_only", "SRM_only", "*"]
    for dv_used in all_dvs_types:
        if dv_used == "*":
            req = "SELECT id, h_p_score, h_a_score FROM solutions_multi_fin"# ORDER BY RANDOM() LIMIT 1000"
            line_thickness = 0.05
            alpha = 0.5
        else:
            req = "SELECT id FROM solutions_multi_fin WHERE h_a > 200e3 AND h_p > 200e3 AND h_a < 500e3 AND h_p < 500e3 AND dv_used = '%s'"%dv_used
            line_thickness = 0.5
            alpha = 1.0
        cur.execute(req)
        res = cur.fetchall()
        ids = [i[0] for i in res]
        # Plot trajectories
        if nice_png_fig:
            plt.figure(figsize=(9,12.73))
        else:
            plt.figure(figsize=(9,5))
        for i, id in enumerate(ids):
            if i % 25 == 0:
                print("Plotting trajectory %i of %i for %s"%(i+1, len(ids), dv_used), end="\r")
            # Load simulation results
            try:
                sim_results = np.load(sys.path[0]+"/optimisation/sim_results/%i.npz"%id)
            except FileNotFoundError:
                continue
            times = sim_results["times"]
            altitudes = sim_results["dep_vars"][:,0]
            if nice_png_fig:
                try:
                    score_a = max(res[i][1], 1/100)
                    score_b = max(res[i][2], 1/100)
                except TypeError:
                    continue
                line_thickness = (min(1/np.sqrt(score_a), 1)/3+min(1/np.sqrt(score_b), 1)*2/3)*0.5
                alpha = (min(1/np.sqrt(score_a), 1)/3+min(1/np.sqrt(score_b), 1)*2/3)*0.75
            plt.plot(times[:1000]/60, altitudes[:1000]/1e3, linewidth=line_thickness, color="C%i"%(i+1), alpha=alpha)
            plt.plot(times[1001:]/60, altitudes[1001:]/1e3, linewidth=line_thickness, color="C%i"%(i+1), alpha=alpha)
            plt.plot(times[998:1002]/60, altitudes[998:1002]/1e3, linewidth=line_thickness, linestyle="solid", color="C%i"%(i+1), alpha=alpha)
        print("Plotted %i trajectories for %s                     "%(len(ids), dv_used))
        # Prettify plot
        plt.xlabel("Time [min]")
        plt.ylabel("Altitude [km]")
        if dv_used == "*":
            if nice_png_fig:
                plt.ylim((-2.5,2750))
                plt.xlim((0,160))
            else:
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
            if nice_png_fig:
                plt.axis('off')
                plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/altitude", transparent=True, dpi=400)
            else:
                plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/altitude.pdf")
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

if pareto_fronts:
    for dv_used in ["init_angle_only", "TVC_only", "SRM_only", "all", "*"]:
        print("Plotting pareto front for %s..."%dv_used)
        # Load entire database into dataframe
        if dv_used == "*":
            fname = "Pareto.pdf"
            req_a = "SELECT h_a_score+h_p_score, * FROM solutions_multi_fin WHERE h_a_score < 15 AND h_p_score < 15"
            req_b = "SELECT * FROM solutions_multi_fin WHERE h_a_score < 25 AND h_p_score < 25 AND mass_score < 25"
        else:
            fname = "Pareto_%s.pdf" % dv_used
            req_a = "SELECT h_a_score+h_p_score, * FROM solutions_multi_fin WHERE h_a_score < 15 AND h_p_score < 15 AND dv_used = '%s' ORDER BY h_a_score+h_p_score+10*mass_score ASC LIMIT 400" % dv_used
            req_b = "SELECT * FROM solutions_multi_fin WHERE h_a_score < 15 AND h_p_score < 15 AND mass_score < 15 AND dv_used = '%s'" % dv_used
        df = pd.read_sql_query(req_a, con)
        plt.figure(figsize=(9, 5))
        # Plot pareto front
        if dv_used in ["init_angle_only", "TVC_only"]:
            fig, ax = plotting.pareto_front(
                x_objective=df["h_a_score"],
                y_objective=df["h_p_score"],
                x_label="Apoapsis altitude score",
                y_label="Periapsis altitude score",
                alpha=0.65
            )
        else:
            fig, ax = plotting.pareto_front(
                x_objective=df["h_a_score"],
                y_objective=df["h_p_score"],
                x_label="Apoapsis altitude score",
                y_label="Periapsis altitude score",
                c_parameter=df["mass_score"],
                c_label="Mass score",
                cmap="viridis",
                alpha=0.65
            )
        # Save the plot
        plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/"+fname)
        # Plot all objectives against each other in 3D
        print("Plotting pareto front in 3D for %s..."%dv_used)
        df = pd.read_sql_query(req_b, con)
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(111, projection="3d")
        r = ax.scatter(df["h_a_score"], df["h_p_score"], df["mass_score"], alpha=0.45, s=1.2, c=df["mass_score"], cmap="jet")
        # Add colormap
        cbar = fig.colorbar(r, ax=ax, shrink=0.6, pad=0.1)
        cbar.set_label("Mass score")
        ax.set_xlabel("Periapsis altitude score")
        ax.set_ylabel("Apoapsis altitude score")
        ax.set_zlabel("Mass score")
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration/"+fname.replace(".", "_3D."))

# Close the database connection
con.close()
