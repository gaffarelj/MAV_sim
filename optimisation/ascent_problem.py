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
import multiprocessing as MP

# Custom imports
from optimisation.run_ascent import MAV_ascent_sim
from thrust.solid_thrust import SRM_thrust
from thrust.models.spherical import spherical_SRM
from thrust.models.multi_fin import multi_fin_SRM

class MAV_problem:

    def __init__(self, design_var_range, thrust_geo_model_1):
        self.design_var_range = design_var_range
        self.thrust_geo_model_1 = thrust_geo_model_1 # (geo model for stage 2 is always spherical)

    def get_bounds(self):
        return self.design_var_range
    
    def get_nobj(self):
        return 2

    def get_nix(self):
        if self.thrust_geo_model_1 == multi_fin_SRM:
            return 1
        else:
            return 0

    def fitness(self):
        # Implement this function because pygmo checks that it exists, even if batch fitness is used instead.
        raise NotImplementedError("This is a batch problem, use batch_fitness instead.")

    def batch_fitness(self, dv_s):
        inputs, fitnesses = [], []

        # Compute number of design variables
        dv_size = len(self.design_var_range[0])
        n_dvs = len(dv_s)//dv_size
        
        # Loop trough the design variables
        for dv_i in range(n_dvs):
            dv = dv_s[dv_i*dv_size:(dv_i+1)*dv_size]

            # Extract data from the design variables:
            launch_angle_1 = dv[0]
            launch_angle_2 = dv[1]
            TVC_angles_y = list(dv[2:7])
            TVC_angles_z = list(dv[7:12])
            SRM_2_geo_vars = dv[12:14]
            SRM_1_geo_vars = dv[14:]

            # Setup SRM model of stage 1
            if self.thrust_geo_model_1 == multi_fin_SRM:
                L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = SRM_1_geo_vars
                R_i_1 = R_i_frac_1 * R_o_1
                L_f = L_f_frac * R_i_1
                w_f = 2*np.pi*(R_i_1-L_f)/N_f*w_f_frac
                SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
            else:
                raise NotImplementedError("Only multi-fin SRM is implemented")

            # Setup SRM model of stage 2
            R_o_2_frac, R_i_2_frac = SRM_2_geo_vars
            R_o_2 = R_o_2_frac * R_o_1
            R_i_2 = R_i_2_frac * R_o_2
            SRM_2_model = spherical_SRM(R_o_2, R_i_2)

            # Set SRM thrust model from geometry models
            SRM_thrust_model_1 = SRM_thrust(SRM_1_model, A_t=0.065, epsilon=45)
            SRM_thrust_model_2 = SRM_thrust(SRM_2_model, A_t=0.005, epsilon=73, p_a=0)

            inputs.append((launch_angle_1, launch_angle_2, TVC_angles_y, TVC_angles_z, SRM_thrust_model_1, SRM_thrust_model_2))

        # Get the fitness by running the simulations in parallel
        with MP.get_context("spawn").Pool(processes=8) as pool:
            outputs = pool.starmap(MAV_ascent_sim, inputs)

        # Return the 1D list of fitnesses
        for output in outputs:
            altitude_score, mass_score = output
            fitnesses.append(altitude_score)
            fitnesses.append(mass_score)

        return fitnesses


if __name__ == "__main__":
    from thrust.models.multi_fin import multi_fin_SRM
    opti_problem = MAV_problem(0, multi_fin_SRM)