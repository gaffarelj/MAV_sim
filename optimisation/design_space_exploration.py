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

# Custom imports
from optimisation.ascent_problem import MAV_problem

# Run only as main  (protect from multiprocessing)
if __name__ == "__main__":
    # Parameters
    n_tot = 20
    n_MC = int(0.35*n_tot)
    n_Sobol = int(0.65*n_tot)

    # Define design variable range
    angle_1_range = [
        [np.deg2rad(30)],
        [np.deg2rad(60)]
    ]
    angle_2_range = [
        [np.deg2rad(60)],
        [np.deg2rad(120)]
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
    for i in range(n_MC):
        for i_dv in range(len(dv_ranges[0])):
            if i_dv in idx_int_des_var:
                dv_s_MC.append(np.random.randint(dv_ranges[0][i_dv], dv_ranges[1][i_dv]))
            else:
                dv_s_MC.append(np.random.uniform(dv_ranges[0][i_dv], dv_ranges[1][i_dv]))

    # [Sobol] Populate design variables following a Sobol sequence
    dv_s_Sobol = []
    if n_Sobol > 0:
        log_n = int(np.ceil(np.log2(n_Sobol)))
        # Take samples from Sobol sequence
        sampler = qmc.Sobol(d=len(dv_ranges[0]), scramble=True, seed=seed)
        samples = sampler.random_base2(m=log_n) # m=3 means 2^3=8 samples between 0 and 1
        # Make values fit in limits
        # x_sobol = x_lims[0] + samples[:, 0] * (x_lims[1] - x_lims[0])
        # y_sobol = y_lims[0] + samples[:, 1] * (y_lims[1] - y_lims[0])
        for i in range(2**log_n):
            for i_dv in range(len(dv_ranges[0])):
                val = dv_ranges[0][i_dv] + samples[i, i_dv] * (dv_ranges[1][i_dv] - dv_ranges[0][i_dv])
                if i_dv in idx_int_des_var:
                    val = int(val)
                dv_s_Sobol.append(val)

    # Combine design variables from MC and Sobol
    dv_s = dv_s_MC + dv_s_Sobol

    # Create problem
    from thrust.models.multi_fin import multi_fin_SRM
    problem = MAV_problem(dv_ranges, multi_fin_SRM)
    print("Running design space exploration with %i inputs (%i MC, %i Sobol)" % (n_MC+n_Sobol, n_MC, n_Sobol))
    fitnesses = problem.batch_fitness(dv_s)
    print(fitnesses)