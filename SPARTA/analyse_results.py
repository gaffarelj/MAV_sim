import numpy as np
import sys
import sys
sys.path.insert(0, "/mnt/c/Users/jerem/OneDrive - Delft University of Technology/current/Thesis/MAV_sim")

# Define the altitudes and satellite names for which there is results to analyse
hs = [100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]

# Define the velocities and densities at the altitude
Vs =   [3503.34,     3490.86,     3478.51,     3466.29,     3454.20,     3442.23,     3430.39,     3418.67,     3407.07,     3384.21,     3361.81,     3339.85,     3318.31]
rhos = [1.82031e-07, 4.29462e-09, 1.94534e-10, 2.15307e-11, 3.60822e-12, 8.64996e-13, 2.74321e-13, 1.01443e-13, 4.26962e-14, 1.10285e-14, 4.70352e-15, 2.97560e-15, 2.29431e-15]

# Define the satellite names and reference areas
sat_names = ["MAV_stage_2"]
ref_areas = [0.144]

# Loop trough the satellite names and the altitudes
for i_s, s_name in enumerate(sat_names):
    for i_h, h in enumerate(hs):
        # Get the sorted file list corresponding to the satellite name and altitude
        try:
            results_data = open(sys.path[0]+"/SPARTA/setup/results_sparta/%s/stats_%skm.dat" % (s_name, h)).readlines()
        # Print a warning and skip this config if no results file could be found
        except FileNotFoundError:
            print("Warning, it seems that the simulation for %s at %skm was not run yet." % (s_name, h))
            continue

        results = []
        reading_results = False
        for res_line in results_data:
            if reading_results:
                results_row = res_line.split()
                if results_row[0] == "Loop":
                    reading_results = False
                else:
                    results.append(results_row)

            if res_line.strip()[:8] == "Step CPU":
                reading_results = True

        results = np.asarray(results, dtype=float)
        try:
            times, fx, fy, fz, ppc = results[:,0], -results[:,-4], -results[:,-3], -results[:,-2], results[:,-1]
        except IndexError:
            print("Warning, it seems that the simulation for %s at %skm was not run properly." % (s_name, h))
            continue
        
        # Print the satellite and altitude
        print("%s at %.1fkm (%i epochs):" % (s_name, h, times[-1]))

        # Print the drag
        print("Drag = %.5e N" % np.mean(fx[-2:]))
        
        # Compute the force coefficients
        cx = 2*np.mean(fx[-3:])/(rhos[i_h]*Vs[i_h]**2*ref_areas[i_s])
        cy = 2*np.mean(fy[-3:])/(rhos[i_h]*Vs[i_h]**2*ref_areas[i_s])
        cz = 2*np.mean(fz[-3:])/(rhos[i_h]*Vs[i_h]**2*ref_areas[i_s])

        # Print the force coefficients
        print("Aero coefficients: [%.5f, %.5f, %.5f]" % (cx, cy, cz))
        if False:
            # Plot the force in each direction
            PU.plot_single(times, fx, "Timestep number [-]", "$F_x$ [N]", "SPARTA/fx_%s_%skm" % (s_name, h))
            PU.plot_single(times, fy, "Timestep number [-]", "$F_y$ [N]", "SPARTA/fy_%s_%skm" % (s_name, h))
            PU.plot_single(times, fz, "Timestep number [-]", "$F_y$ [N]", "SPARTA/fz_%s_%skm" % (s_name, h))

            # Plot the number of particles over time
            PU.plot_single(times, ppc, "Timestep number [-]", "Mean number of particles per cell [-]", "SPARTA/npart_%s_%skm" % (s_name, h))
