import numpy as np
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
import shutil
import tudatpy.util as TU

# Define altitudes
hs = [100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]

tot_epochs = [1500] * len(hs)   # Number of simulation epochs for each altitude (should be multiple of 1000)
run_fractions = [20/30, 10/30]   # Epochs at which to switch from initial run [0] to refinements [1 to -2] to final refinement and run [-1]
refinement_factors = [2]              # Factor by which to scale the grid
refine_region = [False]         # When True, only refine the grid in the region where the vehicle is
# Scale the number of particles by these
particles_scales = {
    100: [1e2, 10],
    125: [5e3, 10],
    150: [5e4, 10],
    175: [1e5, 10],
    200: [5e5, 10],
    225: [1e6, 10],
    250: [5e6, 10],
    275: [1e7, 10],
    300: [5e7, 10],
    350: [1e8, 10],
    400: [5e8, 10],
    450: [1e9, 10],
    500: [5e9, 10]
}
# Set whether to use the exhaust plume or not
use_exhaust = False
# List of vehicle names
vehicle_names = ["MAV_stage_2"]
# List of vehicle reference lengths
L_s = [0.77]
# List of vehicle lengths
L_vehicles = [0.77]
# List of vehicle (semi-)widths to limit the refinement region
widths = [0.5]
# vehicles for which to run what altitude
vehicle_run_h = {h: "MAV_stage_2" for h in hs}

# Define conditions at different orbital altitudes
rhos = [1.82031e-07, 4.29462e-09, 1.94534e-10, 2.15307e-11, 3.60822e-12, 8.64996e-13, 2.74321e-13, 1.01443e-13, 4.26962e-14, 1.10285e-14, 4.70352e-15, 2.97560e-15, 2.29431e-15]
ps =   [4.58579e-03, 1.07168e-04, 7.27184e-06, 9.83263e-07, 2.05882e-07, 6.32226e-08, 2.59061e-08, 1.29177e-08, 7.73174e-09, 4.31813e-09, 3.26268e-09, 2.72788e-09, 2.35800e-09]
Ts =   [128.643,     133.166,     172.404,     178.405,     179.196,     179.330,     179.369,     179.380,     179.382,     179.383,     179.383,     179.382,     179.383]
Vs =   [3503.34,     3490.86,     3478.51,     3466.29,     3454.20,     3442.23,     3430.39,     3418.67,     3407.07,     3384.21,     3361.81,     3339.85,     3318.31]
fracs = [
    np.array([93.270, 2.445, 2.513, 0.965, 0.607, 0.200])/100,
    np.array([81.674, 6.398, 4.972, 3.851, 2.509, 0.596])/100,
    np.array([65.826, 11.923, 4.910, 7.796, 8.533, 1.012])/100,
    np.array([44.654, 16.782, 3.754, 10.801, 22.822, 1.186])/100,
    np.array([24.723, 17.401, 2.297, 11.148, 43.375, 1.055])/100,
    np.array([11.796, 14.477, 1.210, 9.397, 62.345, 0.775])/100,
    np.array([5.102, 10.644, 0.582, 7.033, 76.132, 0.507])/100,
    np.array([1.900, 7.190, 0.249, 4.862, 85.498, 0.301])/100,
    np.array([0.648, 4.581, 0.098, 3.174, 91.333, 0.166])/100,
    np.array([0.074, 1.758, 0.015, 1.271, 96.835, 0.047])/100,
    np.array([0.012, 0.679, 0.003, 0.501, 98.792, 0.014])/100,
    np.array([0.004, 0.280, 0.001, 0.203, 99.507, 0.005])/100,
    np.array([0.003, 0.132, 0.001, 0.088, 99.775, 0.002])/100
]

run_all_cmd = "#!/bin/sh\n"
paraview_surf = ""
paraview_grid = ""
for j, s_name in enumerate(vehicle_names):
    print("\n\n* vehicle", s_name)
    # Create folder for results of this vehicle
    try:
        os.mkdir(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
    except (FileExistsError, OSError):
        try:
            # Un/comment the two following lines to always remove the previous results when new input files are made
            if False:
                shutil.rmtree(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
                os.mkdir(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
        except (PermissionError, OSError):
            print("Warning: could not delete folder", s_name)
    convert_STL = True
    # Loop trough conditions
    for i, h in enumerate(hs):
        particles_scale = particles_scales[h]
        print("\n - With conditions at altitude of %i km:" % h)
        print("     Velocity is of %.2f m/s, and the species are mixed as follows:" % Vs[i])
        print("    ", fracs[i])
        # Inputs
        rho = rhos[i]   # density [kg/m3]
        p = ps[i]       # pressure [Pa]
        T = Ts[i]       # temperature [K]
        u_s = Vs[i]     # free-stream velocity [m/s]
        L = L_s[j]      # reference length [m] (vehicle width)
        h_box = 1.25     # box height [m]
        w_box = 1.25     # box width [m]
        l_box = 1.75    # box length [m]
        # Fraction of each species
        species_frac = fracs[i]
        if round(sum(species_frac), 4) != 1:
            print("Warning, the sum of the species fraction does not add up to 1 (but to %.4f)..." % sum(species_frac)), input("Press ENTER to continue...")

        # Constants
        # Species name, mass, diameter, frontal area
        species_names = ["CO2", "N2", "Ar", "CO", "O", "O2"]
        species_m = np.array([7.31E-26, 4.65E-26, 6.63E-26, 4.65E-26, 2.65E-26, 5.31E-26])
        species_d = np.array([33e-11, 364e-12, 340e-12, 376e-12, 7.4e-11, 346e-12])
        species_sigma = np.pi*species_d**2
        k_b = 1.38e-23  # Boltzmann constant

        # Compute values
        weighted_m = sum(species_frac * species_m)                      # weighted mass of species
        weighted_sigma = sum(species_frac * species_sigma)              # weighted frontal area of species
        nrho = rho / weighted_m                                         # number density [#/m3]
        lambda_f = 1 / (np.sqrt(2) * weighted_sigma * nrho)             # mean free path [m]
        Kn = lambda_f / L                                               # Knudsen number [-]
        # T_ps = weighted_m * u_s**2 / (7*k_b)                            # post-shock T (Cv=2.5*R) [K]
        # cr_ps = np.sqrt(16*k_b*T_ps / (np.pi*weighted_m))               # post-shock average relative velocity [m/s]
        # nrho_ps = 23/3*nrho                                             # post-shock number density [#/m3]
        # lambda_ps = 1 / (np.sqrt(2) * weighted_sigma * nrho_ps)         # post-shock mean free path [m]
        T_ps, cr_ps, nrho_ps, lambda_ps = 4500, 2250, 3e19, 0.05        # Values taken from preliminary SPARTA results
        nu_ps = weighted_sigma * nrho_ps * cr_ps                        # post-shock collision frequency [Hz]
        tau_ps = 1 / nu_ps                                              # post-shock collision time [s]
        dt_mfp = tau_ps / 5                                             # time step [s] (based on mean free path)
        grid_f_mfp = lambda_f / 5                                       # grid dimension before shock [m] (based on mean free path)
        grid_ps_mfp = lambda_ps / 5                                     # post-shock grid dimension [m] (based on mean free path)
        grid_f_vel = u_s*dt_mfp                                         # grid dimension before shock [m] (based on velocity)
        grid_ps_vel = cr_ps*dt_mfp                                      # post-shock grid dimension [m] (based on velocity)
        grid_f = max(min(grid_f_mfp, grid_f_vel, L/50), L/250)          # Take minimum on (or L_ref/50, to avoid grid too small, L_ref/250 to avoid initial grid too big)
        grid_ps = max(min(grid_ps_mfp, grid_ps_vel, L/50), L/250)       # Take minimum grid dimension (or L_ref/50, to avoid grid too small, L_ref/250 to avoid initial grid too big)
        n_real = (nrho + nrho_ps) / 2 * h_box * l_box * w_box           # real number of particles
        n_x = int(l_box / ((grid_f + grid_ps)/2))                       # number of grid segments along x
        n_y = int(w_box / ((grid_f + grid_ps)/2))                       # number of grid segments along y
        n_z = int(h_box / ((grid_f + grid_ps)/2))                       # number of grid segments along z
        n_cells = n_x * n_y * n_z                                       # number of cells
        n_sim = particles_scale[0] * n_cells                            # number of simulated particles (int factor results from analysis to have 10 ppc)
        f_num = n_real / n_sim                                          # f_num for SPARTA

        # Compute the Mach number
        N = 6.0221409e+23
        molar_mass = 8.314 / (weighted_m * N)
        M = u_s / np.sqrt(1.3 * molar_mass * T)
        print("     Mach number equal to %.4f" % M)
        
        # Check that dt is small enough given dx and v
        dx = min(l_box/n_x, w_box/n_y, h_box/n_z)
        dt = min(dt_mfp, dx/u_s*0.75, dx/cr_ps*0.75)                    # Take smallest dt of all (factor of 0.75 to make sure to be below the limit imposed by velocity)
        
        # Compute the accommodation coefficient based on the adsorption of atomic oxygen
        # https://doi.org/10.2514/1.49330
        K = 7.5E-17                     # model fitting parameter
        n_0 = nrho * species_frac[-2]   # number density of the atomic oxygen
        P = n_0 * T                     # atomic oxygen partial pressure
        alpha = K*P/(1+K*P)             # accommodation coefficient
        # test_accommodation = False
        # if test_accommodation and s_name == vehicle_names[0] and h == hs[0]:
        #     P_s = np.linspace(1.5e17, 9e18, 200)
        #     alpha_s = [K*_P/(1+K*_P) for _P in P_s]
        #     PU.plot_single(P_s, alpha_s, "$n_O \cdot T \ [K / m^3]$", "Accommodation coefficient $\\alpha$ [-]", "test_accommodation")

        # Print the results
        print("     Knudsen number is %.5e" % Kn)

        # Only run if specified
        if s_name not in vehicle_run_h[h]:
            print(" - SPARTA input file not created for %s at %ikm, as specified." % (s_name, h))
            continue

        ## Save the results to an input input
        # Convert STL only once
        if convert_STL:
            # Write command to convert surface to ParaView
            paraview_surf += "pvpython ../../tools/surf2paraview.py ../../setup/data/data.%s %s \n" % (s_name, s_name)
            print(" - Converting binary STL to SPARTA surface...")
            with TU.redirect_std():
                # Convert STL from binary to asci
                os.system("python2 \"%s/SPARTA/tools/stl_B2A.py\" \"%s/SPARTA/setup/STL/%s.stl\" -rs" % (sys.path[0], sys.path[0], s_name))
                # Convert STL to data surface for SPARTA
                os.system("python2 \"%s/SPARTA/tools/stl2surf.py\" \"%s/SPARTA/setup/STL/%s_ASCII.stl\" \"%s/SPARTA/setup/data/data.%s\"" % (sys.path[0], sys.path[0], s_name, sys.path[0], s_name))
            convert_STL = False
        print(" - Saving input to file...")

        # Setup the SPARTA inputs
        input_s =  "# SPARTA input file for vehicle %s, for an altitude of %.1fkm\n" % (s_name, h)
        input_s += "print \"\"\nprint \"***** Running SPARTA simulation for %s, at h=%ikm *****\"\nprint \"\"\n" % (s_name, h)
        input_s += "seed                12345\n"
        input_s += "dimension           3\n"
        grid_def = "dimension           3\n"
        input_s += "\n"
        input_s += "global              gridcut 0.1 comm/sort yes surfmax 10000 splitmax 1000\n"
        input_s += "\n"
        input_s += "boundary            o o o\n"
        input_s += "create_box          -%.4f %.4f -%.4f %.4f -%.4f %.4f\n" % (l_box/2, l_box/2, w_box/2, w_box/2, h_box/2, h_box/2)
        grid_def+= "create_box          -%.4f %.4f -%.4f %.4f -%.4f %.4f\n" % (l_box/2, l_box/2, w_box/2, w_box/2, h_box/2, h_box/2)
        input_s += "\n"
        input_s += "create_grid         %i %i %i\n" % (np.ceil(n_x), np.ceil(n_y), np.ceil(n_z))
        input_s += "\n"
        input_s += "balance_grid        rcb part\n"
        input_s += "\n"

        input_s += "global              nrho %.4e fnum %.4e vstream -%.4f 0.0 0.0 temp %.4f\n" % (nrho, f_num, u_s, T)
        input_s += "\n"
        input_s += "species             ../atmo.species CO2 N2 Ar CO O O2\n"
        for n, sp_n in enumerate(species_names):
            input_s += "mixture             atmo %s frac %.4f\n" % (sp_n, species_frac[n])
        if use_exhaust:
            input_s += "mixture             exhaust CO nrho 9.8415e+20 vstream -2370.0 0.0 0.0 temp 2200.0\n" # nrho 9.8415e+22
        input_s += "collide             vss all ../atmo.vss\n"
        input_s += "\n"
        vehicle_centre = L_vehicles[j]/2
        input_s += "read_surf           ../data/data.%s trans %.4f 0 0\n" % (s_name, vehicle_centre)
        input_s += "surf_collide        1 diffuse 293.15 %.4f\n" % (alpha)
        input_s += "surf_modify         all collide 1\n"
        input_s += "\n"
        input_s += "region              vehicle block %.4f %.4f -%.4f %.4f -%.4f %.4f\n" % \
            (-vehicle_centre-0.001, -vehicle_centre+0.001, widths[j], widths[j], widths[j], widths[j])
        input_s += "\n"
        input_s += "fix                 in emit/face atmo xhi zhi zlo yhi ylo\n"
        if use_exhaust:
            input_s += "fix                 prop_out emit/surf exhaust all\n"
        input_s += "\n"
        input_s += "timestep            %.4e\n" % dt
        input_s += "\n"
        input_s += "compute             forces surf all all fx fy fz\n"
        stats_freq = tot_epochs[i]*min(run_fractions)/4
        input_s += "fix                 avg ave/surf all %i %i %i c_forces[*] ave running\n" % (1, stats_freq*2/3, stats_freq)
        input_s += "compute             sum_force reduce sum f_avg[*]\n"
        input_s += "\n"

        # Grid data to save
        grid_data = ["n", "nrho", "massrho", "u"]
        for g_d in grid_data:
            input_s += "compute             %s grid all all %s\n" % (g_d, g_d)
            input_s += "fix                 %s_avg ave/grid all %i %i %i c_%s[*]\n" % (g_d, 1, stats_freq*2/3, stats_freq, g_d)
            input_s += "\n"
        input_s += "compute             avg_ppc reduce ave f_n_avg\n"
        input_s += "\n"
        input_s += "compute             T thermal/grid all all temp\n"
        input_s += "fix                 T_avg ave/grid all %i %i %i c_T[*]\n" % (1, stats_freq*2/3, stats_freq)
        input_s += "\n"
        input_s += "compute             knudsen lambda/grid f_nrho_avg f_T_avg CO2 kall\n"
        input_s += "\n"
        input_s += "stats               %i\n" % stats_freq
        input_s += "stats_style         step cpu wall np nscoll nexit c_sum_force[*] c_avg_ppc\n"
        input_s += "\n"
        input_s += "dump                0 grid all %i ../results_sparta/%s/vals_%ikm_0.*.dat id %s f_T_avg c_knudsen[*]\n" \
            % (stats_freq*2, s_name, h,  " ".join(["f_%s_avg" % _n for _n in grid_data]))
        input_s += "write_grid          ../results_sparta/%s/grid_%ikm_0.dat\n" % (s_name, h)

        grid_def = [grid_def]*len(run_fractions)
        grid_def[0] += "read_grid           ../../setup/results_sparta/%s/grid_%skm_0.dat\n" % (s_name, h)
        # Make sure the run number is a multiple of the stats (to be compatible with the compute/fix)
        run_n = tot_epochs[i] * run_fractions[0]
        if run_n % stats_freq != 0:
            run_n = run_n + (stats_freq - run_n%stats_freq)
        input_s += "run                 %i\n" % run_n
        input_s += "\n"

        new_dt = min(l_box/n_x, w_box/n_y, h_box/n_z)/cr_ps
        for i_refine, epoch_frac in enumerate(run_fractions[1:]):
            input_s += "print \"Refinement level %s\"\n" % (i_refine+1)
            new_dt /= refinement_factors[i_refine]
            # For the new dt, make sure the following condition is satisfied: u_ps*dt < dx
            input_s += "timestep            %.4e\n" % new_dt
            # Refine the grid where the grid Knudsen number is below 5, only in front of the vehicle (coarsen it back when it is above 50)
            input_s += "adapt_grid          all refine coarsen value c_knudsen[2] 5 50 combine min thresh less more cells %i %i %i" \
                 % tuple([refinement_factors[i_refine]]*3)
            # If specified, only refine in the region where the vehicle is
            input_s += " region vehicle one\n" if refine_region[i_refine] else "\n"
            # Increase the number of particles so that the PPC stay > ~10
            input_s += "scale_particles     all %i\n" % particles_scale[i_refine+1]
            f_num /= particles_scale[i_refine+1]
            input_s += "global              fnum %.4e\n" % f_num
            # Make dumps for the new grid
            input_s += "undump              %i\n" % i_refine
            input_s += "dump                %i grid all %i ../results_sparta/%s/vals_%ikm_%i.*.dat id %s f_T_avg c_knudsen[*]\n" \
                % (i_refine+1, stats_freq, s_name, h, i_refine+1, " ".join(["f_%s_avg" % _n for _n in grid_data]))
            # Balance the grid bewteen the processors
            input_s += "balance_grid        rcb part\n"
            input_s += "\n"
            input_s += "write_grid          ../results_sparta/%s/grid_%ikm_%i.dat\n" % (s_name, h, i_refine+1)
            grid_def[i_refine+1] += "read_grid           ../../setup/results_sparta/%s/grid_%skm_%i.dat\n" % (s_name, h, i_refine+1)
            # Make sure the run number is a multiple of the stats (to be compatible with the compute/fix)
            run_n = tot_epochs[i] * epoch_frac
            if run_n % stats_freq != 0:
                run_n = run_n + (stats_freq - run_n%stats_freq)
            input_s += "run                 %i\n" % run_n
            input_s += "\n"
        
        run_all_cmd += "mpirun -np 20 spa_ < in.%s_%skm | tee ../results_sparta/%s/stats_%ikm.dat\n" % (s_name, h, s_name, h)
        for i_r in range(len(run_fractions)):
            paraview_grid += "\n"
            paraview_grid += "rm -rf vals_%s_%skm_%i \n" % (s_name, h, i_r)
            paraview_grid += "rm -rf vals_%s_%skm_%i.pvd \n" % (s_name, h, i_r)
            paraview_grid += "echo 'Converting results of %s at %skm (refinement %i) to ParaView...'\n" % (s_name, h, i_r)
            paraview_grid += "pvpython ../../tools/grid2paraview_original.py def/grid.%s_%skm_%i vals_%s_%skm_%i -r ../../setup/results_sparta/%s/vals_%skm_%i.*.dat \n" % \
                (s_name, h, i_r, s_name, h, i_r, s_name, h, i_r)

        # Write SPARTA inputs to input
        with open(sys.path[0] + "/SPARTA/setup/inputs/in.%s_%skm" % (s_name, h), "w") as input_f:
            input_f.write(input_s)

        # Write grid definition in ParaView folder
        for i_g in range(len(run_fractions)):
            with open(sys.path[0] + "/SPARTA/paraview/grid/def/grid.%s_%skm_%i" % (s_name, h, i_g), "w") as input_f:
                input_f.write(grid_def[i_g])


# Write command to run all SPARTA input files
try:
    with open(sys.path[0] + "/SPARTA/setup/inputs/run_all.sh", "r+") as run_f:
        run_f.seek(0)
        run_f.write(run_all_cmd)
        run_f.truncate()
except FileNotFoundError:
    with open(sys.path[0] + "/SPARTA/setup/inputs/run_all.sh", "w") as run_f:
        run_f.write(run_all_cmd)
# Write command to create ParaView files
paraview_cmd = "#!/bin/sh\n"
paraview_cmd += "cd surf\n"
paraview_cmd += "rm -rf *\n"
paraview_cmd += paraview_surf
paraview_cmd += "cd ../grid\n"
paraview_cmd += paraview_grid
try:
    with open(sys.path[0] + "/SPARTA/paraview/paraview_convert.sh", "r+") as run_f:
        run_f.seek(0)
        run_f.write(paraview_cmd)
        run_f.truncate()
except FileNotFoundError:
    with open(sys.path[0] + "/SPARTA/paraview/paraview_convert.sh", "w") as run_f:
        run_f.write(paraview_cmd)
