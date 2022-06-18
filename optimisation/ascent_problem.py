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
import matplotlib.pyplot as plt
import sqlite3

# Custom imports
from optimisation.run_ascent import MAV_ascent_sim
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust
from thrust.models.spherical import spherical_SRM
from thrust.models.multi_fin import multi_fin_SRM

class MAV_problem:

    def __init__(self, design_var_range, thrust_geo_model_1):
        self.design_var_range = design_var_range
        self.thrust_geo_model_1 = thrust_geo_model_1 # (geo model for stage 2 is always spherical)

    def get_bounds(self):
        return self.design_var_range
    
    def get_nobj(self):
        return 3

    def get_nix(self):
        if self.thrust_geo_model_1 == multi_fin_SRM:
            return 1
        else:
            return 0

    def fitness(self, dv):
        # Return the fitness for a single decision vector
        return self.batch_fitness(dv)

    def batch_fitness(self, dv_s, save_to_db=None):
        inputs, fitnesses = [], []

        # Compute number of design variables
        dv_size = len(self.design_var_range[0])
        n_dvs = len(dv_s)//dv_size
        
        if save_to_db is not None:
            # Connect to the database
            con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db")
            cur = con.cursor()
        
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

            # print("launch_angle_1", np.rad2deg(launch_angle_1))
            # print("launch_angle_2", np.rad2deg(launch_angle_2))
            # print("TVC_angles_y", np.rad2deg(TVC_angles_y))
            # print("TVC_angles_z", np.rad2deg(TVC_angles_z))
            # print("SRM_2_geo_vars", SRM_2_geo_vars)
            # print("SRM_1_geo_vars", SRM_1_geo_vars)

            # Setup SRM model of stage 1
            if self.thrust_geo_model_1 == multi_fin_SRM:
                L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = SRM_1_geo_vars
                R_i_1 = R_i_frac_1 * R_o_1
                L_f = L_f_frac * R_i_1
                w_f = w_f_frac*2*np.pi*(R_i_1-L_f)/N_f
                # Make sure that fins will burn before outer tube
                if w_f/2 > R_o_1 - R_i_1:
                    w_f = 2*(R_o_1 - R_i_1)-1e-7
                SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
            else:
                raise NotImplementedError("Only multi-fin SRM is implemented")

            # Setup SRM model of stage 2
            R_o_2_frac, R_i_2_frac = SRM_2_geo_vars
            R_o_2 = R_o_2_frac * R_o_1
            R_i_2 = R_i_2_frac * R_o_2
            SRM_2_model = spherical_SRM(R_o_2, R_i_2)

            db_id = None
            if save_to_db is not None:
                # Skip if inputs already in database
                req = "SELECT id, h_p_score FROM solutions_multi_fin WHERE "
                req += " AND ".join(["angle_%i = ?"%i for i in range(1,3)])
                req += " AND "
                req += " AND ".join(["TVC_y_%i = ?"%i for i in range(1,6)])
                req += " AND "
                req += " AND ".join(["TVC_z_%i = ?"%i for i in range(1,6)])
                req += " AND "
                req += " AND ".join(["spherical_motor_%i = ?"%i for i in range(1,3)])
                req += " AND "
                req += " AND ".join(["multi_fin_motor_%i = ?"%i for i in range(1,7)])
                cur.execute(req, dv)
                res = cur.fetchall()
                # If results exist...
                if len(res) > 0:
                    db_h_p_score = res[0][1]
                    db_id = res[0][0]
                    # Check if score was already computed...
                    if db_h_p_score is not None:
                        print("Skipping already existing solution with id", db_id)
                        continue
                # Otherwise, add inputs to database
                else:
                    # Set inputs id as latest one + 1
                    highest_id = cur.execute("SELECT MAX(id) FROM solutions_multi_fin").fetchone()[0]
                    if highest_id is None:
                        db_id = 1
                    else:
                        db_id = highest_id + 1
                    # Save the design variables to the database
                    req = "INSERT INTO solutions_multi_fin (id, "
                    req += ", ".join(["angle_%i"%i for i in range(1,3)])
                    req += ", "
                    req += ", ".join(["TVC_y_%i"%i for i in range(1,6)])
                    req += ", "
                    req += ", ".join(["TVC_z_%i"%i for i in range(1,6)])
                    req += ", "
                    req += ", ".join(["spherical_motor_%i"%i for i in range(1,3)])
                    req += ", "
                    req += ", ".join(["multi_fin_motor_%i"%i for i in range(1,7)])
                    req += ", dv_used) VALUES ("
                    req += ", ".join(["?" for i in range(22)])
                    req += ")"
                    cur.execute(req, (db_id,)+tuple(dv)+(save_to_db,))

            # Set SRM thrust model from geometry models
            SRM_thrust_model_1 = SRM_thrust(SRM_1_model, A_t=0.065, epsilon=45)
            SRM_thrust_model_2 = SRM_thrust(SRM_2_model, A_t=0.005, epsilon=73, p_a=0)

            inputs.append((launch_angle_1, launch_angle_2, TVC_angles_y, TVC_angles_z, SRM_thrust_model_1, SRM_thrust_model_2, False, db_id))

        if save_to_db is not None:
            con.commit()
            con.close()

        # Get the fitness by running the simulations in parallel
        with MP.get_context("spawn").Pool(processes=50) as pool:#MP.cpu_count()//2) as pool:
            outputs = pool.starmap(MAV_ascent_sim, inputs)

        # Return the 1D list of fitnesses
        for output in outputs:
            h_p_score, h_a_score, mass_score = output
            fitnesses.append(h_p_score)
            fitnesses.append(h_a_score)
            fitnesses.append(mass_score)

        return fitnesses


if __name__ == "__main__":
    from thrust.models.multi_fin import multi_fin_SRM
    opti_problem = MAV_problem(0, multi_fin_SRM)