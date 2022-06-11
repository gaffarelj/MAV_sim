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
from scipy.stats import qmc
import sqlite3

# Custom imports
from optimisation.ascent_problem import MAV_problem

# Run only as main  (protect from multiprocessing)
if __name__ == "__main__":
    # Parameters
    dv_to_use = "SRM_only" # "all", "init_angle_only", "TVC_only", "SRM_only"
    n_tot = int(14000/0.7)
    n_MC = int(0.3*n_tot)
    n_Sobol = int(0.7*n_tot)
    reset_table = False

    # Define default values for the design variables
    default_angle_1 = [np.deg2rad(57.5)]
    default_angle_2 = [np.deg2rad(90)]
    default_TVC_y = [0, 0.05, 0.1, 0, 0.05]
    default_TVC_z = [0, -0.05, 0.0, 0.05, 0.05]
    default_spherical_geo = [0.165/0.24, 0.0915/0.165]
    default_multi_fin_geo = [1.05, 0.24, 0.175/0.24, 0.05/0.175, 0.02/(2*np.pi*(0.175-0.05)/20), 20]
    default_dv = np.concatenate([default_angle_1, default_angle_2, default_TVC_y, default_TVC_z, default_spherical_geo, default_multi_fin_geo])

    if reset_table:
        # Connect to the database
        con = sqlite3.connect(sys.path[0]+"/optimisation/space_exploration.db")
        cur = con.cursor()
        # Delete the table
        cur.execute("DROP TABLE IF EXISTS solutions_multi_fin")
        # Create the table
        req = "CREATE TABLE solutions_multi_fin (id REAL, angle_1 REAL, angle_2 REAL, "
        req += ", ".join(["TVC_y_%i REAL"%i for i in range(1, 6)])
        req += ", "
        req += ", ".join(["TVC_z_%i REAL"%i for i in range(1, 6)])
        req += ", "
        req += ", ".join(["spherical_motor_%i REAL"%i for i in range(1, 3)])
        req += ", "
        req += ", ".join(["multi_fin_motor_%i REAL"%i for i in range(1, 7)])
        req += ", h_p_score REAL, h_a_score REAL, mass_score REAL, final_time REAL, h_a REAL, h_p REAL, mass REAL, inclination REAL, t_b_1 REAL, t_b_2 REAL, dv_used TEXT)"
        cur.execute(req)

    # Define design variable range
    angle_1_range = [
        [np.deg2rad(30)],   # -> [np.deg2rad(47.5)]
        [np.deg2rad(60)]
    ]
    angle_2_range = [
        [np.deg2rad(60)],   # -> [np.deg2rad(70)]
        [np.deg2rad(120)]   # -> [np.deg2rad(110)]
    ]
    TVC_y_range = [
        [np.deg2rad(-5)]*5,
        [np.deg2rad(5)]*5
    ]
    TVC_z_range = [
        [np.deg2rad(-5)]*5,
        [np.deg2rad(5)]*5
    ]
    spherical_motor_range = [
        [0.3, 0.2],
        [1.0, 0.9]
    ]
    multi_fin_motor_range = [
        [0.3, 0.1, 0.2, 0.25, 0.35, 3],
        [1.25, 0.285, 0.9, 0.75, 0.9, 15]
    ]
    idx_int_des_var = [19]

    # Combine the ranges for all design variables
    dv_ranges = [[], []]
    for dv_range in [angle_1_range, angle_2_range, TVC_y_range, TVC_z_range, spherical_motor_range, multi_fin_motor_range]:
        dv_ranges[0].extend(dv_range[0])
        dv_ranges[1].extend(dv_range[1])

    seed = 123987
    np.random.seed(seed)
    # [MC] Populate design variables at random
    dv_s_MC = []
    n_dvs = len(dv_ranges[0])
    i_start = 0
    if dv_to_use == "init_angle_only":
        n_dvs = 2
        i_start = 0
    elif dv_to_use == "TVC_only":
        n_dvs = 10
        i_start = 2
    elif dv_to_use == "SRM_only":
        n_dvs = 8
        i_start = 12
    for i in range(n_MC):
        dv = default_dv.copy()
        for i_dv in range(n_dvs):
            if i_dv+i_start in idx_int_des_var:
                dv[i_dv+i_start] = np.random.randint(dv_ranges[0][i_dv+i_start], dv_ranges[1][i_dv+i_start])
            else:
                dv[i_dv+i_start] = np.random.uniform(dv_ranges[0][i_dv+i_start], dv_ranges[1][i_dv+i_start])
        dv_s_MC.extend(dv)

    # [Sobol] Populate design variables following a Sobol sequence
    dv_s_Sobol = []
    if n_Sobol > 0:
        log_n = int(np.log2(n_Sobol))
        n_Sobol = 2**log_n
        # Take samples from Sobol sequence
        sampler = qmc.Sobol(d=len(dv_ranges[0]), scramble=True, seed=seed)
        samples = sampler.random_base2(m=log_n) # m=3 means 2^3=8 samples between 0 and 1
        # Make values fit in limits
        # x_sobol = x_lims[0] + samples[:, 0] * (x_lims[1] - x_lims[0])
        # y_sobol = y_lims[0] + samples[:, 1] * (y_lims[1] - y_lims[0])
        for i in range(2**log_n):
            dv = default_dv.copy()
            for i_dv in range(n_dvs):
                val = dv_ranges[0][i_dv+i_start] + samples[i, i_dv+i_start] * (dv_ranges[1][i_dv+i_start] - dv_ranges[0][i_dv+i_start])
                if i_dv+i_start in idx_int_des_var:
                    val = int(val)
                dv[i_dv+i_start] = val
            dv_s_Sobol.extend(dv)

    # Combine design variables from MC and Sobol
    dv_s = dv_s_MC + dv_s_Sobol

    # Create problem
    from thrust.models.multi_fin import multi_fin_SRM
    problem = MAV_problem(dv_ranges, multi_fin_SRM)
    print("Running design space exploration with %i inputs (%i MC, %i Sobol)..." % (n_MC+n_Sobol, n_MC, n_Sobol))
    fitnesses = problem.batch_fitness(dv_s, plot_results=True, save_to_db=dv_to_use)
    print(fitnesses)