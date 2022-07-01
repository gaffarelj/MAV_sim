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
from scipy.stats import qmc

# Tudapty imports
from tudatpy import util

# Custom imports
from optimisation.ascent_problem import MAV_problem

make_plot = False
run_refinement = False      # refinement_ = vary original value at random within 0.1% | refinement2_ = vary original value by 0.1% a random value in its range; 64 samples
refine_mass_opti = False    # refinement3_ = vary original value by 0.25% a random value in its range; 1024 samples

if __name__ == "__main__":

    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
    cur = con.cursor()

    # Get all the optimums (from all optimisations / design space exploration / tuning)
    req = "SELECT * FROM solutions_multi_fin WHERE h_a_score + h_p_score < 0.7 AND mass_score < 1.0 AND dv_used NOT LIKE 'refinement%' AND angle_1 IS NOT NULL and h_a_score IS NOT NULL"
    df = pd.read_sql_query(req, con)

    angle_1_range = [[np.deg2rad(47.5)],[np.deg2rad(60)]]
    angle_2_range = [[np.deg2rad(70)],[np.deg2rad(110)]]
    TVC_z_range = [[np.deg2rad(-5)]*5,[np.deg2rad(5)]*5]
    spherical_motor_range = [[0.3, 0.2],[1.0, 0.9]]
    multi_fin_motor_range = [[0.3, 0.1, 0.2, 0.25, 0.35, 3],[1.25, 0.285, 0.9, 0.75, 0.9, 15]]
    # Combine the ranges for all design variables
    dv_ranges = [[], []]
    for dv_range in [angle_1_range, angle_2_range, TVC_z_range, spherical_motor_range, multi_fin_motor_range]:
        dv_ranges[0].extend(dv_range[0])
        dv_ranges[1].extend(dv_range[1])

    if make_plot:
        plt.figure(figsize=(9,5))

        # Loop trough:
        # - design space exploration
        # - optimiser tuning
        # - optimisation
        # - refinement
        # Plot dots for each (each some color)
        # Plot Pareto front for each (same color as dots, always dotted)
        # Save figures separately and together

        req_dicts = {
            "DSE": ("dv_used = 'all' OR dv_used = 'init_angle_only' OR dv_used = 'TVC_only' OR dv_used = 'SRM_only'", 0.75),
            "Optimiser tuning": ("dv_used LIKE '%_tuning_%'", 0.75),
            "Optimisation": ("dv_used LIKE 'opti_%'", 0.5),
            "Refinement 1": ("dv_used LIKE 'refinement_%' AND dv_used NOT LIKE 'refinement2_%' AND dv_used NOT LIKE 'refinement3_%'", 0.15),
            "Refinement 2": ("dv_used LIKE 'refinement2_%'", 0.1),
            "Refinement 3": ("dv_used LIKE 'refinement3_%'", 0.15)
        }

        fig_all = plt.figure(figsize=(9,5))

        x_fronts, y_fronts = [], []

        for i_dv, dv_used, (dv_req, alpha) in zip(range(len(req_dicts)), req_dicts.keys(), req_dicts.values()):
            print("Plotting for:", dv_used)
            fig_single = plt.figure(figsize=(9,5))
            # Get data for this dv_used
            req = "SELECT * FROM solutions_multi_fin WHERE h_a_score + h_p_score < 0.7 AND mass_score < 1.0 AND h_a_score IS NOT NULL and angle_1 IS NOT NULL AND (%s)"%dv_req
            df_dv = pd.read_sql_query(req, con)
            # Plot dots
            x_objective = df_dv["mass_score"]
            y_objective = df_dv["h_p_score"]+df_dv["h_a_score"]
            optimum_mask = util.pareto_optimums(np.array([x_objective, y_objective]).T)
            x_fronts.extend(x_objective[optimum_mask]), y_fronts.extend(y_objective[optimum_mask])
            # Plot just for this one
            plt.figure(fig_single)
            plt.scatter(x_objective, y_objective, s=10, c="C%i"%i_dv, alpha=0.5)
            plt.step(sorted(x_objective[optimum_mask], reverse=True), sorted(y_objective[optimum_mask], reverse=False), c="C%i"%i_dv, linestyle="dashed", linewidth=1.5)
            plt.xlabel("Mass score [-]")
            plt.ylabel("Altitude score ($h_a + h_p$) [-]")
            plt.grid()
            plt.tight_layout()
            plt.savefig(sys.path[0]+"/plots/optimisation/results/Pareto_%s.pdf"%dv_used.lower().replace(" ", "_"))
            # Plot with all
            plt.figure(fig_all)
            plt.scatter(x_objective, y_objective, s=10, c="C%i"%i_dv, alpha=alpha)
            plt.step(sorted(x_objective[optimum_mask], reverse=True), sorted(y_objective[optimum_mask], reverse=False), c="C%i"%i_dv, linestyle="dashed", label=dv_used, linewidth=1.5)
        # Save plot with all
        x_fronts, y_fronts = np.array(x_fronts), np.array(y_fronts)
        optimum_mask = util.pareto_optimums(np.array([x_fronts, y_fronts]).T)
        plt.step(sorted(x_fronts[optimum_mask], reverse=True), sorted(y_fronts[optimum_mask], reverse=False), c="k", linestyle="dotted", label="From all results", linewidth=1.5)
        plt.xlabel("Mass score [-]")
        plt.ylabel("Altitude score ($h_a + h_p$) [-]")
        plt.grid()
        plt.legend(loc="upper right")
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/results/Pareto_all.pdf")
        plt.close()

    elif run_refinement:

        # Get Pareto front
        points = [
            [df["h_a_score"].iloc[i]+df["h_p_score"].iloc[i], df["mass_score"].iloc[i]]
            for i in range(len(df))
        ]
        optimums = util.pareto_optimums(points)
        print("%i/%i points are Pareto optimums" % (np.count_nonzero(optimums), len(points)))
        optimums_idx = np.where(optimums)[0]
        df = df.iloc[optimums_idx].reset_index()

        # Vary design variable of each optimum following a Sobol sequence within 0.1% of the original values
        variation_bound = 0.1/100
        log_n = 6
        n_Sobol = 2**log_n # 2^6 = 64 samples
        n_dv = 15
        seed = 42

        # Take samples from Sobol sequence
        sampler = qmc.Sobol(d=n_dv, scramble=True, seed=seed)
        samples = sampler.random_base2(m=log_n)
        idx_int_des_var = [14]
        for i_row in range(len(df)):
            dv_s_Sobol = []
            print("*** Running refinement for optimum %i/%i ***" % (i_row+1, len(df)))
            id = df["id"].loc[i_row]
            default_dv = [df["angle_1"].iloc[i_row], df["angle_2"].iloc[i_row], df["TVC_z_1"].iloc[i_row], df["TVC_z_2"].iloc[i_row], df["TVC_z_3"].iloc[i_row], df["TVC_z_4"].iloc[i_row], df["TVC_z_5"].iloc[i_row], df["spherical_motor_1"].iloc[i_row], df["spherical_motor_2"].iloc[i_row], df["multi_fin_motor_1"].iloc[i_row], df["multi_fin_motor_2"].iloc[i_row], df["multi_fin_motor_3"].iloc[i_row], df["multi_fin_motor_4"].iloc[i_row], df["multi_fin_motor_5"].iloc[i_row], df["multi_fin_motor_6"].iloc[i_row]]
            # Make values fit in limits
            for i_sobol in range(2**log_n):
                dv = default_dv.copy()
                for i_dv in range(n_dv):
                    val = dv[i_dv] + variation_bound * (dv_ranges[0][i_dv] + samples[i_sobol, i_dv] * (dv_ranges[1][i_dv] - dv_ranges[0][i_dv]))
                    # if i_dv in idx_int_des_var and not np.isnan(val):
                    #     val = int(val)
                    val = max(dv_ranges[0][i_dv], val)
                    val = min(dv_ranges[1][i_dv], val)
                    dv[i_dv] = val
                if not np.any(np.isnan(dv)):
                    dv_s_Sobol.extend(dv)

            # Create problem
            from thrust.models.multi_fin import multi_fin_SRM
            problem = MAV_problem([[1]*15]*2, multi_fin_SRM)
            fitnesses = problem.batch_fitness(dv_s_Sobol, save_to_db="refinement2_%i"%id)

    elif refine_mass_opti:
        # Get indexes where h_a_score+h_p_score = 0
        idx_0 = np.where(df["h_a_score"]+df["h_p_score"] == 0)[0]
        # Get index where mass_score is min
        idx_mass_opti = np.where(df["mass_score"] == df["mass_score"].iloc[idx_0].min())[0][0]
        id = df["id"].loc[idx_mass_opti]

        # Vary design variable of each optimum following a Sobol sequence within 0.25% of the original values
        variation_bound = 0.25/100
        log_n = 10
        n_Sobol = 2**log_n # 2^10 = 1024 samples
        n_dv = 15
        seed = 42

        # Take samples from Sobol sequence
        sampler = qmc.Sobol(d=n_dv, scramble=True, seed=seed)
        samples = sampler.random_base2(m=log_n)
        idx_int_des_var = [14]
        default_dv = [df["angle_1"].iloc[idx_mass_opti], df["angle_2"].iloc[idx_mass_opti], df["TVC_z_1"].iloc[idx_mass_opti], df["TVC_z_2"].iloc[idx_mass_opti], df["TVC_z_3"].iloc[idx_mass_opti], df["TVC_z_4"].iloc[idx_mass_opti], df["TVC_z_5"].iloc[idx_mass_opti], df["spherical_motor_1"].iloc[idx_mass_opti], df["spherical_motor_2"].iloc[idx_mass_opti], df["multi_fin_motor_1"].iloc[idx_mass_opti], df["multi_fin_motor_2"].iloc[idx_mass_opti], df["multi_fin_motor_3"].iloc[idx_mass_opti], df["multi_fin_motor_4"].iloc[idx_mass_opti], df["multi_fin_motor_5"].iloc[idx_mass_opti], df["multi_fin_motor_6"].iloc[idx_mass_opti]]
        dv_s_Sobol = []
        # Make values fit in limits
        for i_sobol in range(2**log_n):
            dv = default_dv.copy()
            for i_dv in range(n_dv):
                val = dv[i_dv] + variation_bound * (dv_ranges[0][i_dv] + samples[i_sobol, i_dv] * (dv_ranges[1][i_dv] - dv_ranges[0][i_dv]))
                # if i_dv in idx_int_des_var and not np.isnan(val):
                #     val = int(val)
                val = max(dv_ranges[0][i_dv], val)
                val = min(dv_ranges[1][i_dv], val)
                dv[i_dv] = val
            if not np.any(np.isnan(dv)):
                dv_s_Sobol.extend(dv)

        # Create problem
        from thrust.models.multi_fin import multi_fin_SRM
        problem = MAV_problem([[1]*15]*2, multi_fin_SRM)
        fitnesses = problem.batch_fitness(dv_s_Sobol, save_to_db="refinement3_%i"%id)
        
    con.close()