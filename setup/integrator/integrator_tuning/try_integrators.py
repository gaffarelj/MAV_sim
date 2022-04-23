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
    
    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/setup/integrator/integrator_tuning/database.db")
    cur = con.cursor()

    # # /!\ RESET the integrator table /!\
    # cur.execute("DROP TABLE IF EXISTS integrator")
    # cur.execute("CREATE TABLE integrator (method TEXT, coefficients TEXT, end_time REAL, f_evals REAL, error_pos REAL, error_vel REAL, tolerance REAL, dt REAL)")
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
    dts = sorted(np.logspace(-1, -2, num=10), reverse=True)
    inputs = []
    for coefficients in fixed_coefficients:
        for dt in dts:
            # Check in the database that this combination was not already run
            res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=? AND dt=?", ("fixed", str(coefficients).split(".")[-1], dt))
            inputs.append(("fixed", coefficients, dt))
    con.close()
    # coefficients = propagation_setup.integrator.rk_4
    # dt=1.3
    # run_ascent(method="fixed", coefficients=coefficients, dt=dt)
    
    with MP.get_context("spawn").Pool(24) as pool:
        outputs = pool.starmap(run_ascent, inputs)