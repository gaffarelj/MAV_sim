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
    from setup.integrator.run_bench_dt import run_all

    # Define maximum time step to use
    dt = 9.95 # (do no use an int to avoid artifacts with perfect numerical values)
    min_dt = 1e-6

    # Get list of timesteps for which simulations have been run
    filenames = os.listdir(sys.path[0]+"/setup/integrator/benchmark_sim_results")
    filenames.remove(".gitkeep")
    list_dts = sorted([name.replace("dt_", "").replace(".npz", "") for name in filenames])
        
    inputs = []
    while dt > min_dt:
        if "%.4e" % dt not in list_dts:
            inputs.append([dt])
        dt = 10**(np.log10(dt) - 0.1)

    with MP.get_context("spawn").Pool(8) as pool:
        outputs = pool.starmap(run_all, inputs)