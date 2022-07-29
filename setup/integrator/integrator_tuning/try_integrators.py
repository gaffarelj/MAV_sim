if __name__ == "__main__":

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
    import os
    import multiprocessing as MP
    import sqlite3

    # Tudatpy imports
    from tudatpy.kernel.numerical_simulation import propagation_setup

    # Custom imports
    from setup.integrator.integrator_tuning.run_specific_integrator import run_ascent

    # Parameters
    method = "fixed" # "variable" or "fixed"
    dts = sorted(np.logspace(-3, -0.5, num=10), reverse=True)
    tolerances = sorted(np.logspace(-18, -2, num=20), reverse=True)
    allowed_errors = [5e3, 5]
    reset_table = False

    if method not in ["variable", "fixed"]:
        raise ValueError("Method must be 'variable' or 'fixed'")
    
    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/setup/integrator/integrator_tuning/database.db")
    cur = con.cursor()

    # /!\ RESET the integrator table /!\
    if reset_table:
        cur.execute("DROP TABLE IF EXISTS integrator")
        cur.execute("CREATE TABLE integrator (method TEXT, coefficients TEXT, end_time REAL, f_evals REAL, error_pos REAL, error_vel REAL, tolerance REAL, dt REAL)")
        con.commit()

    # # Delete all row where error_pos is -1
    # cur.execute("DELETE FROM integrator WHERE error_pos = -1")
    # cur.execute("DELETE FROM integrator WHERE method = 'variable'")
    # con.commit()

    # Define fixed integrators to try
    fixed_coefficients = [
        propagation_setup.integrator.euler_forward,
        propagation_setup.integrator.rk_4,
        propagation_setup.integrator.explicit_mid_point,
        propagation_setup.integrator.explicit_trapezoid_rule,
        propagation_setup.integrator.ralston,
        propagation_setup.integrator.rk_3,
        propagation_setup.integrator.ralston_3,
        propagation_setup.integrator.SSPRK3,
        propagation_setup.integrator.ralston_4,
        propagation_setup.integrator.three_eight_rule_rk_4
    ]

    # Define variable integrators to try
    variable_coefficients = [
        propagation_setup.integrator.heun_euler,
        propagation_setup.integrator.rkf_12,
        propagation_setup.integrator.rkf_45,
        propagation_setup.integrator.rkf_56,
        propagation_setup.integrator.rkf_78,
        propagation_setup.integrator.rkdp_87,
        propagation_setup.integrator.rkf_89,
        propagation_setup.integrator.rkv_89,
        propagation_setup.integrator.rkf_108,
        propagation_setup.integrator.rkf_1210,
        propagation_setup.integrator.rkf_1412
    ]

    iterate_coefficients = variable_coefficients if method == "variable" else fixed_coefficients#+variable_coefficients
    iterate_input = dts if method == "fixed" else tolerances
    db_col = "dt" if method == "fixed" else "tolerance"

    inputs = []
    N_workers = 44

    #while True:
    for i_inp, inp in enumerate(iterate_input):
        if inp >= 1e-10:
            for i_coeff, coefficients in enumerate(iterate_coefficients):
                res1 = cur.execute("SELECT error_pos, error_vel FROM integrator WHERE method=? AND coefficients=? AND error_pos<? AND error_vel<? AND %s>=?" % db_col, (method, str(coefficients).split(".")[-1], allowed_errors[0]/10, allowed_errors[1]/10, inp)).fetchall()
                res2 = cur.execute("SELECT error_pos, error_vel FROM integrator WHERE method=? AND coefficients=? AND f_evals>=5e4", (method, str(coefficients).split(".")[-1])).fetchall()
                if len(res1) != 0 or len(res2) != 0:
                    pass
                else:
                    # Check in the database that this combination was not already run
                    res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=? AND %s=?" % db_col, (method, str(coefficients).split(".")[-1], inp))
                    if res.fetchone() is None:
                        inputs.append((method, coefficients, inp))
    
    con.close()
    print(inputs)
    input("Press enter to run these inputs...")
    print("Running {} inputs...".format(len(inputs)))
    
    with MP.get_context("spawn").Pool(N_workers) as pool:
        outputs = pool.starmap(run_ascent, inputs)
