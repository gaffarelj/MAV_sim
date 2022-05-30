if __name__ == "__main__":

    import sys
    # Set path to uppermost project level
    sys.path = [p for p in sys.path if p != ""]
    while sys.path[0].split("/")[-1] != "MAV_sim":
        sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

    # Standard imports
    import numpy as np
    import os
    import multiprocessing as MP

    # Custom imports
    from setup.integrator.benchmark.run_bench_dt import run_all

    # Define maximum time step to use
    dt = 9.95e3 # (do no use an int to avoid artifacts with perfect numerical values)
    min_dt = 1e-8
    current_stage = 1
    powered = True
    only_thrust = False
    thrust_rk4 = True

    # Get list of timesteps for which simulations have been run
    filenames = os.listdir(sys.path[0]+"/setup/integrator/benchmark_sim_results")
    list_dts = []
    for name in filenames:
        try:
            if only_thrust:
                list_dts.append(name.replace("thrust_%s_%i_dt_"%("rk4" if thrust_rk4 else "euler", current_stage), "").replace(".npz", ""))
            else:
                list_dts.append(name.replace("%i_%s_dt_"%(current_stage, "V" if powered else "X"), "").replace(".npz", ""))
        except ValueError:
            pass
    list_dts = sorted(list_dts)
        
    inputs = []
    while dt > min_dt:
        if "%.4e" % dt not in list_dts:
            inputs.append([dt, current_stage, powered, only_thrust, True, thrust_rk4])
        dt = 10**(np.log10(dt) - 0.1)

    # Add one more input half the last dt, to compute error in the last dt
    half_dt = inputs[-1][0]*0.49
    if "%.4e" % half_dt not in list_dts:
        inputs.append([half_dt, current_stage, powered, only_thrust, True, thrust_rk4])

    print("Press ENTER to run the following inputs:\n", inputs), input()

    with MP.get_context("spawn").Pool(6) as pool:
        outputs = pool.starmap(run_all, inputs)
