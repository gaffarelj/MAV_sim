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
from scipy.optimize import fsolve
import re

# Tudatpy imports
from tudatpy import plotting
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup

# Custom imports
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.anchor import anchor_SRM
from thrust.models.tubular import tubular_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust
from setup.ascent_framework import MAV_ascent, FakeAeroGuidance

# Connect to the database
con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
cur = con.cursor()

# Parameters
make_plots = False
make_stats = False
plot_best = True
SRM_type = "tubular"

replace_cols = None
if SRM_type == "multi_fin":
    SRM_class = multi_fin_SRM
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "multi_fin_motor_1": "multi_fin_motor_L", "multi_fin_motor_2": "multi_fin_motor_R_o", "multi_fin_motor_3": "multi_fin_motor_R_i", "multi_fin_motor_4": "multi_fin_motor_L_f", "multi_fin_motor_5": "multi_fin_motor_w_f", "multi_fin_motor_6": "multi_fin_motor_N_f"}
    des_vars_names_more = ["multi_fin_motor_L", "multi_fin_motor_R_o", "multi_fin_motor_R_i", "multi_fin_motor_L_f", "multi_fin_motor_w_f", "multi_fin_motor_N_f"]
    des_vars_legends_more = ("$L$", "$R_o$", "$R_i$", "$L_f$", "$w_f$", "$N_f$")
    SRM_range_to_use = [[0.3, 0.1, 0.2, 0.25, 0.35, 3], [1.25, 0.285, 0.9, 0.75, 0.9, 20]]
elif SRM_type == "tubular":
    SRM_class = tubular_SRM
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "tubular_motor_1": "tubular_motor_L", "tubular_motor_2": "tubular_motor_R_o", "tubular_motor_3": "tubular_motor_R_i"}
    des_vars_names_more = ["tubular_motor_L", "tubular_motor_R_o", "tubular_motor_R_i"]
    des_vars_legends_more = ("$L$", "$R_o$", "$R_i$")
    SRM_range_to_use = [[0.3, 0.1, 0.2], [1.25, 0.285, 0.9]]
elif SRM_type == "rod_and_tube":
    SRM_class = rod_and_tube_SRM
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "rod_and_tube_motor_1": "rod_and_tube_motor_L", "rod_and_tube_motor_2": "rod_and_tube_motor_R_o", "rod_and_tube_motor_3": "rod_and_tube_motor_R_mid", "rod_and_tube_motor_4": "rod_and_tube_motor_R_i"}
    des_vars_names_more = ["rod_and_tube_motor_L", "rod_and_tube_motor_R_o", "rod_and_tube_motor_R_mid", "rod_and_tube_motor_R_i"]
    des_vars_legends_more = ("$L$", "$R_o$", "$R_{mid}$", "$R_i$")
    SRM_range_to_use = [[0.3, 0.1, 0.2, 0.2], [1.25, 0.285, 0.9, 0.9]]
elif SRM_type == "anchor":
    SRM_class = anchor_SRM
    replace_cols = {"dv_used": "generation", "spherical_motor_1": "spherical_motor_R_o", "spherical_motor_2": "spherical_motor_R_i", "anchor_motor_1": "anchor_motor_L", "anchor_motor_2": "anchor_motor_R_o", "anchor_motor_3": "anchor_motor_R_i", "anchor_motor_4": "anchor_motor_N_a", "anchor_motor_5": "anchor_motor_w", "anchor_motor_6": "anchor_motor_r_f", "anchor_motor_7": "anchor_motor_delta_s"}
    des_vars_names_more = ["anchor_motor_L", "anchor_motor_R_o", "anchor_motor_R_i", "anchor_motor_N_a", "anchor_motor_w", "anchor_motor_r_f", "anchor_motor_delta_s"]
    des_vars_legends_more = ("$L$", "$R_o$", "$R_i$", "$w$", "$r_f$", "$\delta_s$", "$N_a$")
    SRM_range_to_use = [[0.3, 0.05, 0.2, 0.3, 0.05, 0.1, 2], [1.25, 0.285, 0.6, 0.85, 0.95, 0.75, 6]]
else:
    raise ValueError("Unknown SRM type")


def get_df_gens(df):
    # Extract generation number from dv_used column, and replace column name
    df["dv_used"] = df["dv_used"].str.split("_").str[-1].astype(int)
    df.rename(columns=replace_cols, inplace=True)
    return df

des_vars_legends = [("$\\theta_1$", "$\\theta_2$"), ("$TVC_{z,1}$", "$TVC_{z,2}$", "$TVC_{z,3}$", "$TVC_{z,4}$", "$TVC_{z,5}$"), ("$R_o$", "$R_i$"), des_vars_legends_more]
des_vars_legends_flat = [item for tup in des_vars_legends for item in tup]
des_vars_names = ["angle_1", "angle_2", "TVC_z_1", "TVC_z_2", "TVC_z_3", "TVC_z_4", "TVC_z_5", "spherical_motor_R_o", "spherical_motor_R_i"] + des_vars_names_more

if make_plots:
    for seed in [42, 13, 123, 846, 579]:
        # Load optimisation results
        req = "SELECT * FROM solutions_%s WHERE dv_used LIKE 'opti_%s_%%' AND h_p_score IS NOT NULL AND h_p_score < 5e2 AND h_a_score < 5e2" % (SRM_type, seed)
        df = pd.read_sql_query(req, con)
        df = get_df_gens(df)

        # Define the range in which the design variables can vary
        launch_angle_1_range = np.deg2rad([47.5, 60])
        launch_angle_2_range = np.deg2rad([70, 110])
        TVC_range = np.deg2rad([-5, 5])
        N_TVC_nodes = 5
        spherical_SRM_range = [[0.3, 0.2], [1.0, 0.9]]
        design_var_range = (
            [launch_angle_1_range[0], launch_angle_2_range[0], *[TVC_range[0]]*N_TVC_nodes, *spherical_SRM_range[0], *SRM_range_to_use[0]],
            [launch_angle_1_range[1], launch_angle_2_range[1], *[TVC_range[1]]*N_TVC_nodes, *spherical_SRM_range[1], *SRM_range_to_use[1]]
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
        plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/history_mm_%i.pdf" % (SRM_type, seed))

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
        plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/Pareto_%i.pdf" % (SRM_type, seed))
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
        plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/history_std_%i.pdf" % (SRM_type, seed))

    ## Get optimum values
    # Select where h_comb_score = 0 and mass_score is min
    req = "SELECT * FROM solutions_%s WHERE dv_used LIKE 'opti_%%' AND h_p_score+h_a_score = 0 ORDER BY mass_score ASC LIMIT 3" % SRM_type
    df_opt = pd.read_sql(req, con=con)
    # Order df_opt by mass_score
    df_opt = df_opt.sort_values(by="mass_score")
    for i_row in range(df_opt.shape[0]):
        # Get the design variables
        print("Initial angles [deg]: %.2f, %.2f" % tuple([np.rad2deg(df_opt.iloc[i_row][col_n]) for col_n in ["angle_1", "angle_2"]]))
        print("TVC angles: %.2f, %.2f, %.2f, %.2f, %.2f" % tuple([np.rad2deg(df_opt.iloc[i_row]["TVC_z_%i"%i]) for i in range(1, 6)]))

        # Setup SRM model of stage 1
        if SRM_class == multi_fin_SRM:
            SRM_name = "multi_fin"
            L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = [df_opt.iloc[i_row]["%s_motor_%i"%(SRM_type, i)] for i in range(1, 7)]
            R_i_1 = R_i_frac_1 * R_o_1
            L_f = L_f_frac * R_i_1
            w_f = w_f_frac*2*np.pi*(R_i_1-L_f)/N_f
            # Make sure that fins will burn before outer tube
            if w_f/2 > R_o_1 - R_i_1:
                w_f = 2*(R_o_1 - R_i_1)-1e-7
            SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
        elif SRM_class == tubular_SRM:
            SRM_name = "tubular"
            L, R_o_1, R_i_frac_1 = [df_opt.iloc[i_row]["%s_motor_%i"%(SRM_type, i)] for i in range(1, 4)]
            R_i_1 = R_i_frac_1 * R_o_1
            SRM_1_model = tubular_SRM(R_o_1, R_i_1, L)
        elif SRM_class == rod_and_tube_SRM:
            SRM_name = "rod_and_tube"
            L, R_o_1, R_mid_frac, R_i_frac_1 = [df_opt.iloc[i_row]["%s_motor_%i"%(SRM_type, i)] for i in range(1, 5)]
            R_mid = R_mid_frac * R_o_1
            R_i_1 = R_i_frac_1 * R_mid
            SRM_1_model = rod_and_tube_SRM(R_o_1, R_mid, R_i_1, L)
        elif SRM_class == anchor_SRM:
            SRM_name = "anchor"
            L, R_o_1, R_i_frac_1, w_frac, r_f_frac, delta_s_frac, N_a = [df_opt.iloc[i_row]["%s_motor_%i"%(SRM_type, i)] for i in range(1, 8)]
            N_a = int(N_a)
            R_i_1 = R_i_frac_1 * R_o_1
            w = w_frac * (R_o_1 - R_i_1) / 3
            r_f = r_f_frac * (R_o_1 - 3 * w - R_i_1) / 2
            delta_s = delta_s_frac * 2 * R_i_1 * np.sin(np.pi/N_a)
            # Make sure that R_i is always valid
            if np.arcsin( (delta_s + 2 * w)/(2 * (R_i_1 + w)) ) + np.arcsin( (r_f + w)/(R_i_1 + 2 * w + r_f) ) >= np.pi/N_a:
                print("Uh oh")
                R_i_1 = fsolve(lambda x: np.arcsin( (delta_s + 2 * w)/(2 * (x + w)) ) + np.arcsin( (r_f + w)/(x + 2 * w + r_f) ) - np.pi/N_a, R_i_1)[0]+1e-5
            SRM_1_model = anchor_SRM(R_o_1, R_i_1, N_a, w, r_f, delta_s, L, run_checks=False)
        else:
            raise NotImplementedError("SRM model not implemented for optimisation")

        R_o_2_frac, R_i_2_frac = [df_opt.iloc[i_row]["spherical_motor_%i"%i] for i in range(1, 3)]
        R_o_2 = R_o_2_frac * R_o_1
        R_i_2 = R_i_2_frac * R_o_2
        SRM_2_model = spherical_SRM(R_o_2, R_i_2)

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

if make_stats:
    req = "SELECT * FROM solutions_%s WHERE angle_1 IS NOT NULL and h_a_score IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY mass_score ASC LIMIT 100"%SRM_type
    df_opt = pd.read_sql(req, con=con)
    # Only keep column with des_vars_names and scores
    cols_to_keep = ["angle_%i"%i for i in range(1, 3)]
    cols_to_keep += ["TVC_z_%i"%i for i in range(1, 6)]
    cols_to_keep += ["spherical_motor_%i"%i for i in range(1, 3)]
    cols_to_keep += ["%s_motor_%i"%(SRM_type, i) for i in range(1, len(SRM_range_to_use[0])+1)]
    cols_to_keep += ["mass_score"]
    df_opt = df_opt[cols_to_keep]
    # Replace column names
    df_opt.rename(columns=replace_cols, inplace=True)
    # Convert first 7 columns to degrees
    df_stats = df_opt.describe()
    print(df_stats.iloc[1:,:7].apply(lambda x: np.rad2deg(x)))
    print(df_stats.iloc[1:,7:])
    # Print std as a percentage of mean
    print(df_stats.iloc[2,:]/df_stats.iloc[1,:]*100)
    # Take median optimum
    df_median = df_stats.iloc[1,:-1]
    print("Median optimum:")
    print(df_median)
    
if plot_best:
    # req = "SELECT * FROM solutions_%s WHERE angle_1 IS NOT NULL and h_a_score IS NOT NULL AND h_a_score + h_p_score = 0 ORDER BY mass_score ASC LIMIT 1" % SRM_type
    req = "SELECT * FROM solutions_%s WHERE h_a_score = 0 AND h_p > 305e3 ORDER BY mass_score ASC LIMIT 1" % SRM_type
    df_best = pd.read_sql(req, con=con).iloc[0]
    print(df_best)
    
    # Setup SRM model of stage 1
    if SRM_class == multi_fin_SRM:
        SRM_name = "multi_fin"
        L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = df_best[["%s_motor_%i"%(SRM_type, i) for i in range(1, 7)]]
        R_i_1 = R_i_frac_1 * R_o_1
        L_f = L_f_frac * R_i_1
        w_f = w_f_frac*2*np.pi*(R_i_1-L_f)/N_f
        # Make sure that fins will burn before outer tube
        if w_f/2 > R_o_1 - R_i_1:
            w_f = 2*(R_o_1 - R_i_1)-1e-7
        SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
    elif SRM_class == tubular_SRM:
        SRM_name = "tubular"
        L, R_o_1, R_i_frac_1 = df_best[["%s_motor_%i"%(SRM_type, i) for i in range(1, 4)]]
        R_i_1 = R_i_frac_1 * R_o_1
        SRM_1_model = tubular_SRM(R_o_1, R_i_1, L)
    elif SRM_class == rod_and_tube_SRM:
        SRM_name = "rod_and_tube"
        L, R_o_1, R_mid_frac, R_i_frac_1 = df_best[["%s_motor_%i"%(SRM_type, i) for i in range(1, 5)]]
        R_mid = R_mid_frac * R_o_1
        R_i_1 = R_i_frac_1 * R_mid
        SRM_1_model = rod_and_tube_SRM(R_o_1, R_mid, R_i_1, L)
    elif SRM_class == anchor_SRM:
        SRM_name = "anchor"
        L, R_o_1, R_i_frac_1, w_frac, r_f_frac, delta_s_frac, N_a = df_best[["%s_motor_%i"%(SRM_type, i) for i in range(1, 8)]]
        N_a = int(N_a)
        R_i_1 = R_i_frac_1 * R_o_1
        w = w_frac * (R_o_1 - R_i_1) / 3
        r_f = r_f_frac * (R_o_1 - 3 * w - R_i_1) / 2
        delta_s = delta_s_frac * 2 * R_i_1 * np.sin(np.pi/N_a)
        # Make sure that R_i is always valid
        if np.arcsin( (delta_s + 2 * w)/(2 * (R_i_1 + w)) ) + np.arcsin( (r_f + w)/(R_i_1 + 2 * w + r_f) ) >= np.pi/N_a:
            print("Uh oh")
            R_i_1 = fsolve(lambda x: np.arcsin( (delta_s + 2 * w)/(2 * (x + w)) ) + np.arcsin( (r_f + w)/(x + 2 * w + r_f) ) - np.pi/N_a, R_i_1)[0]+1e-5
        SRM_1_model = anchor_SRM(R_o_1, R_i_1, N_a, w, r_f, delta_s, L, run_checks=False)
    else:
        raise NotImplementedError("SRM model not implemented for optimisation")


    SRM_thrust_model_1 = SRM_thrust(SRM_1_model, A_t=0.065, epsilon=45)
    R_o_2_frac, R_i_2_frac = df_best[["spherical_motor_%i"%i for i in range(1, 3)]]
    R_o_2 = R_o_2_frac * R_o_1
    R_i_2 = R_i_2_frac * R_o_2
    SRM_stage_2 = spherical_SRM(R_o_2, R_i_2)
    SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)
    mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
    mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p
    body_fixed_thrust_direction_z = [
        [df_best["TVC_z_%i"%i] for i in range(1, 6)],
        0
    ]

    # print("Angles:", np.rad2deg([df_best["angle_%i"%i] for i in range(1, 3)]))
    # print("TVC:", np.rad2deg([df_best["TVC_z_%i"%i] for i in range(1, 6)]))
    # print("Spherical SRM DV:", R_o_2_frac, R_i_2_frac)
    # print("Spherical SRM:", R_o_2, R_i_2)
    # print("Anchor SRM DV:", L, R_o_1, R_i_frac_1, w_frac, r_f_frac, delta_s_frac, N_a)
    # print("Anchor SRM:", R_o_1, R_i_1, N_a, w, r_f, delta_s, L)
    # print("Rod and tube SRM DV:", L, R_o_1, R_mid_frac, R_i_frac_1)
    # print("Rod and tube SRM:", R_o_1, R_mid, R_i_1, L)
    print("Mass:", mass_1, mass_2)
    # exit()

    ascent_model = MAV_ascent(
        launch_epoch = 0,
        launch_lat = np.deg2rad(18.85),
        launch_lon = np.deg2rad(77.52),
        launch_h = -2.5e3,
        mass_stages = [mass_1, mass_2],            
        launch_angles = [df_best["angle_%i"%i] for i in range(1, 3)],
        thrust_models = [SRM_thrust_model_1, SRM_thrust_model_2],
        target_orbit_h = 300e3,
        target_orbit_i = np.deg2rad(25),
        max_a = 15 * 9.80665,
        max_AoA = np.deg2rad(4),
        body_fixed_thrust_direction_y=[0, 0],
        body_fixed_thrust_direction_z=body_fixed_thrust_direction_z
    )
    stage_res = []
    t_b_1, t_b_2 = 0, 0
    for stage in [1, 2]:
        ascent_model.create_bodies(stage=stage, add_sun=True, use_new_coeffs=True, custom_exponential_model=True)
        ascent_model.create_accelerations(use_cpp=True, better_precision=True)
        guidance_object = FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, ascent_model.current_body, silence_warnings=True)
        ascent_model.create_initial_state()
        ascent_model.create_dependent_variables_to_save()
        ascent_model.dependent_variables_to_save.append(propagation_setup.dependent_variable.geodetic_latitude(ascent_model.current_name, "Mars"))
        ascent_model.dependent_variables_to_save.append(propagation_setup.dependent_variable.longitude(ascent_model.current_name, "Mars"))
        ascent_model.create_termination_settings(end_time=160*60, cpu_time_termination=30)
        ascent_model.create_propagator_settings()
        ascent_model.create_integrator_settings(better_accuracy=False)
        times, states, dep_vars = ascent_model.run_simulation()
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,1])
        if stage == 1:
            t_b_1 = ascent_model.thrust.burn_time
            t_sep = times[-1]
            if final_h < 0:
                break
        else:
            t_b_2 = ascent_model.thrust.burn_time

        burn_times = ascent_model.thrust_models[stage-1].saved_burn_times
        magnitudes = ascent_model.thrust_models[stage-1].saved_magnitudes
        ind_below_0 = [i for i in range(len(magnitudes)) if magnitudes[i] < 0]
        idx_0 = min(ind_below_0) if len(ind_below_0) != 0 else -1
        burn_times = burn_times[:idx_0]
        magnitudes = magnitudes[:idx_0]
        burn_times.append(burn_times[-1])
        magnitudes.append(0)
        # plt.figure(figsize=(4.5,3.25))
        # plt.plot(burn_times, np.asarray(magnitudes)/1e3)
        # plt.xlabel("Burn time [s]"), plt.ylabel("Thrust [kN]")
        # plt.grid(), plt.tight_layout()
        # plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/SRM%i_thrust.pdf"%(SRM_type,stage))
        # plt.close()
        print("Burn time:", burn_times[-1])
        
    if stage == 1:
        # Extract results from first propagation only if stage 2 was not used
        times = stage_res[0][0]
        states = stage_res[0][1]
        dep_vars = stage_res[0][2]
    else:
        # Combine results from both propagations
        times = np.concatenate((stage_res[0][0], stage_res[1][0]))
        states = np.concatenate((stage_res[0][1], stage_res[1][1]))
        dep_vars = np.concatenate((stage_res[0][2], stage_res[1][2]))
    times = (times - times[0])/60
    flight_path_angles = np.rad2deg(dep_vars[:,0])
    altitudes = dep_vars[:,1]
    force_coeffs = dep_vars[:,2:5]
    airspeeds = dep_vars[:,5]
    tot_accs = dep_vars[:,6]
    mach_numbers = dep_vars[:,7]
    mass = dep_vars[:,8]
    angle_of_attacks = np.rad2deg(dep_vars[:,9])
    positions = dep_vars[:,10:13]
    a_SH = dep_vars[:,13]
    a_thrust = dep_vars[:,14]
    a_aero = dep_vars[:,15]
    dyna_pressures = dep_vars[:,16]
    velocities = dep_vars[:,17:20]
    full_a_thrust = dep_vars[:,20:23]
    try:
        idx_crop = np.where(times >= t_sep/60+3)[0][0]
        idx_crop_zoom = np.where(times >= 9.5)[0][0]
        idx_crop_no_margin = np.where(times >= t_sep/60+t_b_2/60)[0][0]
    except IndexError:
        idx_crop = -1

    latitudes = dep_vars[:,-2]
    longitudes = dep_vars[:,-1]
    R = 3389.5
    ground_distance = R * (np.sin(longitudes - longitudes[0]) + np.cos(latitudes - latitudes[0]))
    ground_distance = ground_distance - ground_distance[0]

    # Plot SRM geometry and thrust profile
    ax = SRM_1_model.plot_geometry(add_title=False)
    plt.tight_layout()
    ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    ax.grid(False)
    ax.set(frame_on=False)
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/SRM1.pdf"%SRM_type)
    plt.close()
    ax = SRM_stage_2.plot_geometry(add_title=False)
    plt.tight_layout()
    ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    ax.grid(False)
    ax.set(frame_on=False)
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/SRM2.pdf"%SRM_type)
    plt.close()

    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 7))
    # ax1.plot(times, altitudes/1e3)
    # ax1.set_xlabel("Time [min]")
    # ax1.set_ylabel("Altitude [km]")
    # ax1.grid()
    # ax2.plot(times, airspeeds)
    # ax2.set_xlabel("Time [min]")
    # ax2.set_ylabel("Airspeed [m/s]")
    # ax2.grid()
    # ax3.plot(ground_distance[:idx_crop_no_margin], altitudes[:idx_crop_no_margin]/1e3)
    # ax3.set_xlabel("Ground distance [km]")
    # ax3.set_ylabel("Altitude [km]")
    # ax3.set_aspect("equal")
    # ax3.grid()
    # ax4.plot(times[:idx_crop], mass[:idx_crop])
    # ax4.set_xlabel("Time [min]")
    # ax4.set_ylabel("Mass [kg]")
    # ax4.grid()
    # plt.tight_layout()
    # plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/optimum.pdf"%SRM_type)

    plt.figure(figsize=(4.5, 3))
    plt.plot(times, altitudes/1e3)
    plt.grid()
    plt.xlabel("Time since launch [min]")
    plt.ylabel("Altitude [km]")
    plt.ylim(-10, 360)
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/altitude.pdf"%SRM_type)

    plt.figure(figsize=(9, 3.5))
    plt.plot(ground_distance[:idx_crop_no_margin], altitudes[:idx_crop_no_margin]/1e3)
    plt.grid()
    plt.xlabel("Ground distance [km]")
    plt.ylabel("Altitude [km]")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/2D.pdf"%SRM_type)

    # fig = plt.figure(figsize=(7, 6))
    # ax = fig.add_subplot(111, projection="3d")
    # ax.scatter(states[:,0], states[:,1], states[:,2])
    # ax.set_xlabel("X [km]")
    # ax.set_ylabel("Y [km]")
    # ax.set_zlabel("Z [km]")
    # plt.tight_layout()
    # plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/3D.pdf"%SRM_type)

    # plt.figure(figsize=(9, 5))
    # plt.plot(times[:idx_crop_zoom], altitudes[:idx_crop_zoom]/1e3)
    # plt.grid()
    # plt.xlabel("Time since launch [min]")
    # plt.ylabel("Altitude [km]")
    # plt.tight_layout()
    # plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/altitude_zoomed.pdf"%SRM_type)

    plt.figure(figsize=(4.5, 3))
    plt.plot(times, airspeeds)
    plt.grid()
    plt.xlabel("Time since launch [min]")
    plt.ylabel("Airspeed [m/s]")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/airspeed.pdf"%SRM_type)

    plt.figure(figsize=(4.5, 3))
    plt.plot(times[:idx_crop], mass[:idx_crop])
    plt.grid()
    plt.xlabel("Time since launch [min]")
    plt.ylabel("Rocket mass [kg]")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/mass.pdf"%SRM_type)

    plt.figure(figsize=(4.5, 3))
    plt.plot(times[:idx_crop], flight_path_angles[:idx_crop])
    plt.grid()
    plt.xlabel("Time since launch [min]")
    plt.ylabel("Flight path angle [deg]")
    plt.tight_layout()
    plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/flight_path.pdf"%SRM_type)

    # plt.figure(figsize=(9, 5))
    # plt.plot(times[:idx_crop], tot_accs[:idx_crop], label="Total", linestyle="dotted", color="black")
    # plt.plot(times[:idx_crop], a_SH[:idx_crop], label="SH D/O 4")
    # plt.plot(times[:idx_crop], a_thrust[:idx_crop], label="Thrust")
    # plt.plot(times[:idx_crop], a_aero[:idx_crop], label="Aerodynamic")
    # plt.grid()
    # plt.xlabel("Time since launch [min]")
    # plt.ylabel("Acceleration [m/s$^2$]")
    # plt.legend()
    # # plt.yscale("symlog", linthresh=1)
    # plt.tight_layout()
    # plt.savefig(sys.path[0]+"/plots/optimisation/optimum/%s/accelerations.pdf"%SRM_type)

con.close()


# TODO: plot group of DV for objectives = 0 for h; < 0.95 for mass, to show a good (and normally robust) design. Make one such plot per set of DV