if __name__ == "__main__":

    import sys
    # Add tudatpy path
    sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
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
    dt = 9.85 # (do no use an int to avoid artifacts with perfect numerical values)

    # Get list of timesteps for which simulations have been run
    filenames = os.listdir(sys.path[0]+"/setup/integrator/benchmark_sim_results")
    filenames.remove(".gitkeep")
    list_dts = sorted([float(name.replace("dt_", "").replace(".npz", "")) for name in filenames])

    if len(list_dts) != 0:
        dt = list_dts[0]
        
    inputs = []
    for i in range(100):
        inputs.append([dt])
        dt = 10**(np.log10(dt) - 0.1)

    with MP.get_context("spawn").Pool(processes=int(MP.cpu_count()//2)) as pool:
        outputs = pool.starmap(run_all, inputs)