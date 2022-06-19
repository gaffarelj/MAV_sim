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
from matplotlib import pyplot as plt

# Connect to the database
con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db")
cur = con.cursor()

algorithms = ["NSGA2", "MACO", "NSPSO"]

# Get all results for each algorithm
for algorithm in algorithms:
    # Get all results for this algorithm
    cur.execute("SELECT h_p_score, h_a_score, mass_score, dv_used FROM solutions_multi_fin WHERE dv_used LIKE '%s_%%_tuning_%%'" % algorithm)
    results = cur.fetchall()
    h_p_scores, h_a_scores, mass_scores, generations = {}, {}, {}, []
    for res in results:
        gen = int(res[3].split("_")[-1])
        if gen not in generations:
            generations.append(gen)
            h_p_scores[gen] = [res[0]]
            h_a_scores[gen] = [res[1]]
            mass_scores[gen] = [res[2]]
        else:
            if res[0] < 1e3 and res[1] < 1e3 and res[2] < 1e3:
                h_p_scores[gen].append(res[0])
                h_a_scores[gen].append(res[1])
                mass_scores[gen].append(res[2])
    generations = np.asarray(generations)

    boxplots = False

    if not boxplots:
        plt.figure(figsize=(9,5))
        # Plot the full history
        h_p_scores_mins, h_a_scores_mins, mass_scores_mins = [], [], []
        for gen in sorted(generations):
            plt.scatter(np.ones(len(h_p_scores[gen]))*(gen-0.25), h_p_scores[gen], s=3.0, alpha=0.45, color="C0", marker="v")
            plt.scatter(np.ones(len(h_a_scores[gen]))*gen, h_a_scores[gen], s=3.0, alpha=0.45, color="C1", marker="^")
            plt.scatter(np.ones(len(mass_scores[gen]))*(gen+0.25), mass_scores[gen], s=3.0, alpha=0.45, color="C2", marker="P")
            h_p_scores_mins.append(min(h_p_scores[gen]))
            h_a_scores_mins.append(min(h_a_scores[gen]))
            mass_scores_mins.append(min(mass_scores[gen]))
        # Plot the mins
        plt.scatter(generations-0.25, h_p_scores_mins, color="C0", s=15.0, label="$h_p$ score", marker="v")
        plt.scatter(generations, h_a_scores_mins, color="C1", s=15.0, label="$h_a$ score", marker="^")
        plt.scatter(generations+0.25, mass_scores_mins, color="C2", s=15.0, label="mass score", marker="P")
        plt.ylabel("Score")
        plt.yscale("symlog", linthresh=0.5)
        for x in np.arange(0.5, max(generations)-0.4, 1):
            plt.axvline(x, color="k", alpha=0.25, linestyle="--", linewidth=0.5)
        plt.grid(axis="y")
        # highest_max = max(max(h_p_scores_mins), max(h_a_scores_mins), max(mass_scores_mins))
        # plt.ylim(-highest_max/75, highest_max*1.25)
        plt.legend()
    else:
        # Plot the results on a 3x1 plot
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
        ax1.boxplot(h_p_scores.values())
        ax2.boxplot(h_a_scores.values())
        ax3.boxplot(mass_scores.values())
        # Set y lims to -0.1 to 100
        ax1.set_ylim(-1, 150)
        ax2.set_ylim(-0.5, 25)
        ax3.set_ylim(-0.1, 2.5)
        # Prettify plot
        ax1.set_ylabel("$h_p$ score [-]"), ax2.set_ylabel("$h_a$ score [-]"), ax3.set_ylabel("Mass score [-]")
        ax1.grid(), ax2.grid(), ax3.grid()

    plt.xlabel("Generation")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/tuning/%s_history.pdf" % algorithm)

    # TODO: make same plot for design variables