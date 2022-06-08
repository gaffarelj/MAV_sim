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
            line_thickness = 0.15
        else:
            req = "SELECT id FROM solutions_multi_fin WHERE h_a > 200e3 AND h_p > 200e3 AND h_a < 500e3 AND h_p < 500e3 AND dv_used = '%s'"%dv_used
            line_thickness = 0.5
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
            plt.plot(times[:1000]/60, altitudes[:1000]/1e3, linewidth=line_thickness, color="C%i"%(i+1))
            plt.plot(times[1001:]/60, altitudes[1001:]/1e3, linewidth=line_thickness, color="C%i"%(i+1))
            plt.plot(times[998:1002]/60, altitudes[998:1002]/1e3, linewidth=line_thickness, linestyle="dashed", color="C%i"%(i+1))
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
            plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration_altitude.pdf")
        else:
            plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration_altitude_%s_sampled.pdf"%dv_used)
        plt.close()

if analyse_correlation:
    for dv_used in ["init_angle_only", "TVC_only", "SRM_only"]:#, "all", "*"]:
        if dv_used == "init_angle_only":
            vars = ["h_a", "h_p", "angle_1", "angle_2"]
        elif dv_used == "TVC_only":
            vars = ["h_a", "h_p"]
            [vars.append("TVC_y_%i"%i) for i in range(1,6)]
            [vars.append("TVC_z_%i"%i) for i in range(1,6)]
        elif dv_used == "SRM_only":
            vars = ["h_a", "h_p", "mass"]
            [vars.append("spherical_motor_%i"%i) for i in range(1,3)]
            [vars.append("multi_fin_motor_%i"%i) for i in range(1,7)]
        else:
            raise NotImplementedError
        # Load entire database into dataframe
        req = "SELECT %s FROM solutions_multi_fin WHERE h_a_score < 1e2 AND h_p_score < 1e2 AND dv_used = '%s'"%(", ".join(vars), dv_used)
        df = pd.read_sql_query(req, con)
        # # Make pairplot
        # print("Making pairplots for %s..."%dv_used)
        # sns.pairplot(df, kind="kde")
        # plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration_pairplot_%s.pdf"%dv_used)#, dpi=200)
        # plt.close()
        # Make heatmap
        print("Making heatmap for %s..."%dv_used)
        # plt.figure(figsize=(25, 25))
        sns.heatmap(df.corr(), annot=True, fmt=".2f", cmap="vlag")
        plt.savefig(sys.path[0]+"/plots/optimisation/design_space_exploration_heatmap_%s.pdf"%dv_used)#, dpi=200)
        plt.tight_layout()
        plt.close()

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