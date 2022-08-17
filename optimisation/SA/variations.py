import sys, os

from setup.integrator.benchmark.run_bench_dt import resample
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
import sqlite3
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import multiprocessing as MP
import seaborn as sns
import pandas as pd
from scipy.interpolate import interp1d

# Tudat imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice
from tudatpy import util

# Custom imports
from setup.ascent_framework import MAV_ascent, FakeAeroGuidance
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust
from thrust.models.anchor import anchor_SRM
from thrust.models.spherical import spherical_SRM

create_table = False
run_error_altitude = False
run_error_latlon = False
run_error_sepdelay = False
run_dens_var = False
try_MCD_winds = False
run_error_SRM = False
run_misalignment = False
run_payload_var = False
plot_all = True

N = 10000

np.random.seed(42)
baseline_h_p, baseline_h_a, baseline_i = 309556.6911757747, 341604.1105327816, 0.3348063892742567

def apo_peri_altitude_from_state(cartesian_state):
    final_keplerian_state = element_conversion.cartesian_to_keplerian(
        cartesian_state,
        spice.get_body_gravitational_parameter("Mars")
    )
    final_a, final_e, final_i = final_keplerian_state[:3]
    R_Mars = spice.get_average_radius("Mars")
    h_peri, h_apo = final_a * (1 - final_e) - R_Mars, final_a * (1 + final_e) - R_Mars
    # print(h_peri, h_apo, final_i)
    return h_peri, h_apo, final_i

def run_sim(id, start_h, angle_1, angle_2, TVC_angles_z, SRM1_geo, SRM2_geo, start_lat=np.deg2rad(18.85), \
            start_lon=np.deg2rad(77.52), sep_delay=0, dens_f=1.0, MCD_winds=None, SRM_size_error=0, misalign=0, payload_mass_var=0):
    R_o_1, R_i_1, N_a, w, r_f, delta_s, L = SRM1_geo
    R_o_2, R_i_2 = SRM2_geo
    R_i_1, w, r_f, delta_s, R_i_2 = R_i_1 + SRM_size_error, w + SRM_size_error, r_f + SRM_size_error, delta_s + SRM_size_error, R_i_2 + SRM_size_error
    SRM_1_model = anchor_SRM(R_o_1, R_i_1, N_a, w, r_f, delta_s, L, run_checks=False)
    SRM_2_model = spherical_SRM(R_o_2, R_i_2)
    thrust_model_1 = SRM_thrust(SRM_1_model, A_t=0.065, epsilon=45)
    thrust_model_2 = SRM_thrust(SRM_2_model, A_t=0.005, epsilon=73, p_a=0)

    # for stage in [1, 2]:
    #     if stage == 1:
    #         times, magnitudes, *_ = thrust_model_1.simulate_full_burn(dt=(2.0e-6,True), use_cpp=True)
    #         mdots = thrust_model_1.saved_m_dot_s
    #     elif stage == 2:
    #         times, magnitudes, *_, masses = thrust_model_2.simulate_full_burn(dt=1.5e-2, use_cpp=True)
    #         mdots = thrust_model_2.saved_m_dot_s
    #     times, magnitudes, mdots = np.asarray(times), np.asarray(magnitudes), np.asarray(mdots)
    #     np.savez(sys.path[0]+"/data/opti_thrust_%i.npz"%stage, times=times, magnitudes=magnitudes, masses=mdots)
    # exit()

    mass_2 = 47.5 + thrust_model_2.M_innert + thrust_model_2.M_p
    mass_1 = 65 + mass_2 + thrust_model_1.M_innert + thrust_model_1.M_p + payload_mass_var

    ascent_model = MAV_ascent(
        launch_epoch = 0,
        launch_lat = start_lat,
        launch_lon = start_lon,
        launch_h = start_h,
        mass_stages = [mass_1, mass_2],            
        launch_angles = [angle_1, angle_2],
        thrust_models = [thrust_model_1, thrust_model_2],
        target_orbit_h = 300e3,
        target_orbit_i = np.deg2rad(25),
        max_a = 15 * 9.80665,
        max_AoA = np.deg2rad(4),
        body_fixed_thrust_direction_y=[ 0, 0],
        body_fixed_thrust_direction_z=[ TVC_angles_z, 0 ]
    )

    # Setup and run simulation for both stages
    stage_res = []
    better_accuracy = False
    for stage in [1, 2]:
        ascent_model.create_bodies(stage=stage, add_sun=True, use_new_coeffs=True, custom_exponential_model=True, dens_f=dens_f, use_MCD_winds_only=MCD_winds)
        if SRM_size_error == 0:
            ascent_model.create_accelerations(thrust_fname=sys.path[0]+"/data/opti_thrust_%i.npz"%stage, sep_delay=sep_delay, misalign=misalign)
        else:
            ascent_model.create_accelerations(use_cpp=(stage==1), better_precision=True)
        guidance_object = FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, ascent_model.current_body, silence_warnings=True)
        ascent_model.create_initial_state()
        ascent_model.create_dependent_variables_to_save(False)
        ascent_model.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(ascent_model.current_name, "Mars"))
        ascent_model.create_termination_settings(end_time=160*60, cpu_time_termination=30, exact_time=True)
        ascent_model.create_propagator_settings()
        ascent_model.create_integrator_settings(better_accuracy=better_accuracy)
        # with util.redirect_std():
        times, states, dep_vars = ascent_model.run_simulation()
        stage_res.append([times, states, dep_vars])
        final_h = max(dep_vars[:,0])
        if stage == 1:
            if final_h < 0:
                break
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
    
    h_peri, h_apo, inclination = apo_peri_altitude_from_state(states[-1,:6])
    h_p_diff, h_a_diff, i_diff = baseline_h_p-h_peri, baseline_h_a-h_apo, baseline_i-inclination
    print("periapsis error: %.2f km, apoapsis error: %.2f km, inclination error: %.2e deg"%(h_p_diff/1e3, h_a_diff/1e3, np.rad2deg(i_diff)))
    
    con = sqlite3.connect(sys.path[0]+"/optimisation/SA/SA.db", timeout=30)
    cur = con.cursor()
    cur.execute("UPDATE initial_variations SET h_p_diff = ?, h_a_diff = ?, i_diff = ? WHERE id = ?", (h_p_diff, h_a_diff, i_diff, id))
    con.commit()
    con.close()

    return h_p_diff, h_a_diff, i_diff

if __name__ == "__main__":
    con = sqlite3.connect(sys.path[0]+"/optimisation/SA/SA.db", timeout=30)
    cur = con.cursor()
    if create_table:
        cur.execute("DROP TABLE IF EXISTS initial_variations")
        req = "CREATE TABLE initial_variations (id REAL, param_used_1 TEXT, param_used_2 TEXT, param_value_1 REAL, param_value_2 REAL, h_p_diff REAL, h_a_diff REAL, i_diff REAL)"
        cur.execute(req)
        con.commit()
        exit()

    if not plot_all:
        R_o_2_frac, R_i_2_frac = 0.757257442766104, 0.538750686956536
        L, R_o_1, R_i_frac_1, w_frac, r_f_frac, delta_s_frac, N_a = 0.981308268066878, 0.243814374516337, 0.498146161472015, 0.763186548669054, 0.454733829237075, 0.510452525008623, 4.83731772484036
        angle_1, angle_2 = 0.886429558218858, 1.55566726316892
        TVC_angles_z = [-0.0109554560643352, -0.0175445546050466, -0.00231919058881665, 0.0329777177194789, 0.030065869849508]

        # Lat x lon: 6.6 x 7.7 km https://mars.nasa.gov/resources/25491/perseverance-rover-landing-ellipse-in-jezero-crater/
        # Lat x lon: 0.11157 x 0.13016 deg
        # Alt: 2500-2600 (figure to put in mission overview)
        launch_lats, launch_lons = [np.deg2rad(18.85)], [np.deg2rad(77.52)]
        altitudes = [-2.55e3]
        sep_delays = [0]
        dens_variations = [1.0]
        MCD_winds = [None]
        SRM_errors = [0]
        thrust_misalignments = [0]
        payload_mass_variations = [0]

        if run_error_altitude:
            altitudes = np.random.normal(-2550, 50, size=N)
        elif run_error_latlon:
            launch_lats = np.random.normal(np.deg2rad(18.85), np.deg2rad(0.11157), size=N)
            launch_lons = np.random.normal(np.deg2rad(77.52), np.deg2rad(0.13016), size=N)
        elif run_error_sepdelay:
            sep_delays = np.random.gamma(1, 2, N)
        elif run_dens_var: # https://www.researchgate.net/publication/223688876_The_influence_of_geomagnetic_and_solar_variabilities_on_lower_thermosphere_density
            dens_variations = np.random.normal(1.0, 0.35, N) # max 1.34
        elif try_MCD_winds:
            MCD_Ds = np.arange(0, 23.9, 1.0)
            MCD_Ls = np.arange(0, 359, 1.0)
            MCD_winds = np.array(np.meshgrid(MCD_Ls, MCD_Ds)).T.reshape(-1, 2)
        elif run_error_SRM:
            SRM_errors = np.random.normal(0, 3e-4, N)
        elif run_misalignment:
            thrust_misalignments = np.random.normal(0, np.deg2rad(0.3), N)
        elif run_payload_var:
            payload_mass_variations = np.random.normal(0, 0.25, N)

        N_a = int(N_a)
        R_i_1 = R_i_frac_1 * R_o_1
        w = w_frac * (R_o_1 - R_i_1) / 3
        r_f = r_f_frac * (R_o_1 - 3 * w - R_i_1) / 2
        delta_s = delta_s_frac * 2 * R_i_1 * np.sin(np.pi/N_a)
        # Make sure that R_i is always valid
        if np.arcsin( (delta_s + 2 * w)/(2 * (R_i_1 + w)) ) + np.arcsin( (r_f + w)/(R_i_1 + 2 * w + r_f) ) >= np.pi/N_a:
            R_i_1 = fsolve(lambda x: np.arcsin( (delta_s + 2 * w)/(2 * (x + w)) ) + np.arcsin( (r_f + w)/(x + 2 * w + r_f) ) - np.pi/N_a, R_i_1)[0]+1e-5
        R_o_2 = R_o_2_frac * R_o_1
        R_i_2 = R_i_2_frac * R_o_2

        SRM1_geo = R_o_1, R_i_1, N_a, w, r_f, delta_s, L
        SRM2_geo = R_o_2, R_i_2
        # run_sim(9876543210, altitudes[0], angle_1, angle_2, TVC_angles_z, SRM1_geo, SRM2_geo, launch_lats[0], launch_lons[0], SRM_size_error=0.001), exit()
            
        plt.figure(figsize=(9,5))
        inputs = []
        for start_h in altitudes:
            for i_ll in range(len(launch_lats)):
                for sep_delay in sep_delays:
                    for dens_f in dens_variations:
                        for MCD_wind in MCD_winds:
                            for SRM_error in SRM_errors:
                                for misalign in thrust_misalignments:
                                    for payload_var in payload_mass_variations:
                                        launch_lat, launch_lon = launch_lats[i_ll], launch_lons[i_ll]
                                        SRM1_geo = R_o_1, R_i_1, N_a, w, r_f, delta_s, L
                                        SRM2_geo = R_o_2, R_i_2
                                        # Check if already run
                                        if run_error_altitude:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("start_h", start_h))
                                        elif run_error_latlon:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 = ? AND param_value_1 = ? AND param_value_2 = ?", ("start_lat", "start_lon", launch_lat, launch_lon))
                                        elif run_error_sepdelay:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("sep_delay", sep_delay))
                                        elif run_dens_var:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("dens_var", dens_f))
                                        elif try_MCD_winds:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 = ? AND param_value_1 = ? AND param_value_2 = ?", ("MCD_D", "MCD_Ls", MCD_wind[0], MCD_wind[1]))
                                        elif run_error_SRM:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("SRM_error", SRM_error))
                                        elif run_misalignment:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("misalign", misalign))
                                        elif run_payload_var:
                                            cur.execute("SELECT id, h_p_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 IS NULL AND param_value_1 = ?", ("payload_var", payload_var))

                                        res = cur.fetchall()
                                        h_p_diff_db = None
                                        if len(res) != 0:
                                            h_p_diff_db = res[0][1]
                                            last_id = res[0][0]-1
                                        else:
                                            # Get last db id
                                            cur.execute("SELECT MAX(id) FROM initial_variations")
                                            last_id = cur.fetchone()[0]
                                            if last_id is None:
                                                last_id = 1
                                            # Insert new values
                                            if run_error_altitude:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "start_h", start_h))
                                            elif run_error_latlon:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_used_2, param_value_1, param_value_2) VALUES (?, ?, ?, ?, ?)", (last_id+1, "start_lat", "start_lon", launch_lat, launch_lon))
                                            elif run_error_sepdelay:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "sep_delay", sep_delay))
                                            elif run_dens_var:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "dens_var", dens_f))
                                            elif try_MCD_winds:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_used_2, param_value_1, param_value_2) VALUES (?, ?, ?, ?, ?)", (last_id+1, "MCD_D", "MCD_Ls", MCD_wind[0], MCD_wind[1]))
                                            elif run_error_SRM:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "SRM_error", SRM_error))
                                            elif run_misalignment:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "misalign", misalign))
                                            elif run_payload_var:
                                                cur.execute("INSERT INTO initial_variations (id, param_used_1, param_value_1) VALUES (?, ?, ?)", (last_id+1, "payload_var", payload_var))
                                        # Add to inputs
                                        if h_p_diff_db is None:
                                            inputs.append((last_id+1, start_h, angle_1, angle_2, TVC_angles_z, SRM1_geo, SRM2_geo, launch_lat, launch_lon, sep_delay, dens_f, MCD_wind, SRM_error, misalign, payload_var))
                                        else:
                                            print("Skipping, already run with ID %i." % last_id)
        con.commit()
        # Run sims in parallel in batch of 500
        while len(inputs) > 0:
            if len(inputs) >= 500:
                inputs_current = inputs[:500]
                inputs = inputs[500:]
            else:
                inputs_current = inputs
                inputs = []
            with MP.get_context("spawn").Pool(MP.cpu_count()//2*0+50) as pool:
                pool.starmap(run_sim, inputs_current)

    else:
        ### Initial altitude error
        cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND i_diff < 0.0007 AND i_diff > -0.0007 AND param_value_1 >= -2675 AND param_value_1 <= -2425", ("start_h",))
        res = cur.fetchall()
        altitudes = np.array(res)[:,0]
        altitudes_diffs = altitudes + 2550
        peri_errors, apo_errors, incli_errors = np.array(res)[:,1]/1e3, np.array(res)[:,2]/1e3, np.rad2deg(np.array(res)[:,3])
        print("Plotting for %i results"%len(altitudes))

        # Make dataframe
        df = pd.DataFrame({"Launch altitude error [m]":altitudes_diffs, "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/altitude_pairplot.pdf")

        ### Initial lat/lon error
        cur.execute("SELECT param_value_1, param_value_2, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND param_used_2 = ? AND h_p_diff IS NOT NULL AND i_diff < 0.0025 AND i_diff > -0.0025 AND param_value_1 IS NOT NULL AND param_value_2 IS NOT NULL", ("start_lat","start_lon",))
        res = cur.fetchall()
        lats = np.array(res)[:,0]
        lat_diffs = np.rad2deg(lats) - 18.85
        lons = np.array(res)[:,1]
        lon_diffs = np.rad2deg(lons) - 77.52
        peri_errors = -np.array(res)[:,2]/1e3
        apo_errors = -np.array(res)[:,3]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,4])
        print("Plotting for %i results"%len(lats))

        # Make dataframe
        df = pd.DataFrame({"Launch latitude error [deg]":lat_diffs, "Launch longitude error [deg]":lon_diffs, "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/latlon_pairplot.pdf")

        ### Stage separation delay
        # cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND i_diff < 0.0025 AND i_diff > -0.0025 AND param_value_1 IS NOT NULL", ("sep_delay",))
        cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND param_value_1 IS NOT NULL", ("sep_delay",))
        res = cur.fetchall()
        sep_delays = np.array(res)[:,0]
        peri_errors = -np.array(res)[:,1]/1e3
        apo_errors = -np.array(res)[:,2]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,3])
        print("Plotting for %i results"%len(sep_delays))

        # Make dataframe
        df = pd.DataFrame({"Stage separation delay [s]":sep_delays, "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/sepdelay_pairplot.pdf")

        ### Density
        cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND param_value_1 IS NOT NULL", ("dens_var",))
        res = cur.fetchall()
        dens_vars = np.array(res)[:,0]
        peri_errors = -np.array(res)[:,1]/1e3
        apo_errors = -np.array(res)[:,2]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,3])
        print("Plotting for %i results"%len(dens_vars))

        # Make dataframe
        df = pd.DataFrame({"Density scaling [-]":dens_vars, "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/densvar_pairplot.pdf")

        ### Winds
        cur.execute("SELECT h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND param_value_1 IS NOT NULL", ("MCD_D",))
        res = cur.fetchall()
        peri_errors = -np.array(res)[:,0]/1e3
        apo_errors = -np.array(res)[:,1]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,2])
        print("Plotting for %i results"%len(peri_errors))

        # Make dataframe
        df = pd.DataFrame({"Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/MCD_pairplot.pdf")

        ### SRM Size
        cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND param_value_1 IS NOT NULL", ("SRM_error",))
        res = cur.fetchall()
        SRM_errors = np.array(res)[:,0]*1e3
        peri_errors = -np.array(res)[:,1]/1e3
        apo_errors = -np.array(res)[:,2]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,3])
        print("Plotting for %i results"%len(SRM_errors))

        # Make dataframe
        df = pd.DataFrame({"SRM sizing error [mm]":SRM_errors, "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/SRM_pairplot.pdf")

        ### SRM Size
        cur.execute("SELECT param_value_1, h_p_diff, h_a_diff, i_diff FROM initial_variations WHERE param_used_1 = ? AND h_p_diff IS NOT NULL AND param_value_1 IS NOT NULL", ("misalign",))
        res = cur.fetchall()
        misalignments = np.array(res)[:,0]
        peri_errors = -np.array(res)[:,1]/1e3
        apo_errors = -np.array(res)[:,2]/1e3
        incli_errors = -np.rad2deg(np.array(res)[:,3])
        print("Plotting for %i results"%len(misalignments))

        # Make dataframe
        df = pd.DataFrame({"Thrust misalignment [deg]":np.rad2deg(misalignments), "Periapsis error [km]":peri_errors, "Apoapsis error [km]":apo_errors, "Inclination error [deg]":incli_errors})
        plt.figure(figsize=(9,9))
        g = sns.PairGrid(df, diag_sharey=False)
        g.map_diag(sns.kdeplot, shade=True)
        g.map_upper(sns.kdeplot, shade=True)
        g.map_lower(sns.scatterplot)
        g.tight_layout()
        plt.savefig(sys.path[0]+"/plots/optimisation/SA/initial/misalign_pairplot.pdf")

    con.close()

# TODO: replot all errors * -1 and adapt discussion (+ figure caption)