import sys, os
from thrust.models.anchor import anchor_SRM

from thrust.models.tubular import tubular_SRM
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
from scipy.optimize import fsolve

# Custom imports
from optimisation.run_ascent import MAV_ascent_sim
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust
from thrust.models.spherical import spherical_SRM
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.tubular import tubular_SRM
from thrust.models.anchor import anchor_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM

# Tudatpy imports
from tudatpy import plotting

class MAV_problem:

    def __init__(self, design_var_range, thrust_geo_model_1, save_to_db=None):
        self.design_var_range = design_var_range
        self.thrust_geo_model_1 = thrust_geo_model_1 # (geo model for stage 2 is always spherical)
        self.save_to_db = save_to_db
        self.i_gen = None

    def get_bounds(self):
        return self.design_var_range
    
    def get_nobj(self):
        return 3

    def get_nix(self):
        if self.thrust_geo_model_1 in [multi_fin_SRM, anchor_SRM]:
            return 1
        else:
            return 0

    def fitness(self, dv):
        # Return the fitness for a single decision vector
        return self.batch_fitness(dv)

    def batch_fitness(self, dv_s, save_to_db=None, plot_geo=False, plot_thrust=False, better_accuracy=False):
        inputs, fitnesses = [], []
        if self.i_gen is None:
            self.i_gen = 0
        else:
            self.i_gen += 1

        if save_to_db is None:
            save_to_db = self.save_to_db

        if save_to_db is not None and "~N~" in save_to_db:
            save_to_db = save_to_db.replace("~N~", str(self.i_gen))

        # Compute number of design variables
        dv_size = len(self.design_var_range[0])
        n_dvs = len(dv_s)//dv_size
        
        # Connect to the database
        con = sqlite3.connect(sys.path[0]+"/optimisation/design_space.db", timeout=30)
        cur = con.cursor()
        
        # Loop trough the design variables
        i_exist, scores_exist = [], []
        for dv_i in range(n_dvs):
            dv = dv_s[dv_i*dv_size:(dv_i+1)*dv_size]

            # Extract data from the design variables:
            launch_angle_1 = dv[0]
            launch_angle_2 = dv[1]
            # TVC_angles_y = list(dv[2:7])
            TVC_angles_z = list(dv[2:7])
            SRM_2_geo_vars = dv[7:9]
            SRM_1_geo_vars = dv[9:]

            # print("launch_angle_1", np.rad2deg(launch_angle_1))
            # print("launch_angle_2", np.rad2deg(launch_angle_2))
            # print("TVC_angles_y", np.rad2deg(TVC_angles_y))
            # print("TVC_angles_z", np.rad2deg(TVC_angles_z))
            # print("SRM_2_geo_vars", SRM_2_geo_vars)
            # print("SRM_1_geo_vars", SRM_1_geo_vars)

            # Setup SRM model of stage 1
            if self.thrust_geo_model_1 == multi_fin_SRM:
                SRM_name = "multi_fin"
                L, R_o_1, R_i_frac_1, L_f_frac, w_f_frac, N_f = SRM_1_geo_vars
                R_i_1 = R_i_frac_1 * R_o_1
                L_f = L_f_frac * R_i_1
                w_f = w_f_frac*2*np.pi*(R_i_1-L_f)/N_f
                # Make sure that fins will burn before outer tube
                if w_f/2 > R_o_1 - R_i_1:
                    w_f = 2*(R_o_1 - R_i_1)-1e-7
                SRM_1_model = multi_fin_SRM(R_o_1, R_i_1, int(N_f), w_f, L_f, L)
            elif self.thrust_geo_model_1 == tubular_SRM:
                SRM_name = "tubular"
                L, R_o_1, R_i_frac_1 = SRM_1_geo_vars
                R_i_1 = R_i_frac_1 * R_o_1
                SRM_1_model = tubular_SRM(R_o_1, R_i_1, L)
            elif self.thrust_geo_model_1 == rod_and_tube_SRM:
                SRM_name = "rod_and_tube"
                L, R_o_1, R_mid_frac, R_i_frac_1 = SRM_1_geo_vars
                R_mid = R_mid_frac * R_o_1
                R_i_1 = R_i_frac_1 * R_mid
                SRM_1_model = rod_and_tube_SRM(R_o_1, R_mid, R_i_1, L)
            elif self.thrust_geo_model_1 == anchor_SRM:
                SRM_name = "anchor"
                L, R_o_1, R_i_frac_1, w_frac, r_f_frac, delta_s_frac, N_a = SRM_1_geo_vars
                N_a = int(N_a)
                R_i_1 = R_i_frac_1 * R_o_1
                w = w_frac * (R_o_1 - R_i_1) / 3
                r_f = r_f_frac * (R_o_1 - 3 * w - R_i_1) / 2
                delta_s = delta_s_frac * 2 * R_i_1 * np.sin(np.pi/N_a)
                # Make sure that R_i is always valid
                if np.arcsin( (delta_s + 2 * w)/(2 * (R_i_1 + w)) ) + np.arcsin( (r_f + w)/(R_i_1 + 2 * w + r_f) ) >= np.pi/N_a:
                    R_i_1 = fsolve(lambda x: np.arcsin( (delta_s + 2 * w)/(2 * (x + w)) ) + np.arcsin( (r_f + w)/(x + 2 * w + r_f) ) - np.pi/N_a, R_i_1)[0]+1e-5
                SRM_1_model = anchor_SRM(R_o_1, R_i_1, N_a, w, r_f, delta_s, L, run_checks=False)
            else:
                raise NotImplementedError("SRM model not implemented for optimisation")

            # Setup SRM model of stage 2
            R_o_2_frac, R_i_2_frac = SRM_2_geo_vars
            R_o_2 = R_o_2_frac * R_o_1
            R_i_2 = R_i_2_frac * R_o_2
            SRM_2_model = spherical_SRM(R_o_2, R_i_2)

            db_id = None
            skip_sim = False
            # Skip if inputs already in database
            req = "SELECT id, h_p_score, h_a_score, mass_score, dv_used FROM solutions_" + SRM_name + " WHERE "
            req += " AND ".join(["angle_%i = ?"%i for i in range(1,3)])
            # req += " AND "
            # req += " AND ".join(["TVC_y_%i = ?"%i for i in range(1,6)])
            req += " AND "
            req += " AND ".join(["TVC_z_%i = ?"%i for i in range(1,6)])
            req += " AND "
            req += " AND ".join(["spherical_motor_%i = ?"%i for i in range(1,3)])
            req += " AND "
            req += " AND ".join(["%s_motor_%i = ?"%(SRM_name, i) for i in range(1,len(SRM_1_geo_vars)+1)])
            cur.execute(req, dv)
            res = cur.fetchall()
            # If results exist...
            if len(res) > 0 and not better_accuracy:
                db_id = res[0][0]
                db_h_p_score = res[0][1]
                db_h_a_score = res[0][2]
                db_mass_score = res[0][3]
                # Check if score was already computed...
                if db_h_p_score is not None:
                    print("Skipping DV %i; solution exists with id %i" % (dv_i, db_id))
                    i_exist.append(dv_i)
                    scores_exist.append([db_h_p_score, db_h_a_score, db_mass_score])

                    if save_to_db is not None:
                        # Check that these scores are not yet associated with a database entry with the same dv_used
                        req = "SELECT id FROM solutions_" + SRM_name + " WHERE h_p_score = ? AND h_a_score = ? AND mass_score = ? AND dv_used = ?"
                        cur.execute(req, [db_h_p_score, db_h_a_score, db_mass_score, save_to_db])
                        res = cur.fetchall()
                        # Save in db as if sim was run
                        if len(res) == 0:
                            highest_id = cur.execute("SELECT MAX(id) FROM solutions_" + SRM_name).fetchone()[0]
                            if highest_id is None:
                                db_id = 1
                            else:
                                db_id = highest_id + 1
                            # Save the design variables to the database
                            req = "INSERT INTO solutions_" + SRM_name + " (id, h_p_score, h_a_score, mass_score, dv_used) VALUES (?, ?, ?, ?, ?)"
                            cur.execute(req, (db_id, db_h_p_score, db_h_a_score, db_mass_score, save_to_db))
                    skip_sim = True
            # Otherwise, add inputs to database
            elif save_to_db is not None and not skip_sim:
                # Set inputs id as latest one + 1
                highest_id = cur.execute("SELECT MAX(id) FROM solutions_" + SRM_name).fetchone()[0]
                if highest_id is None:
                    db_id = 1
                else:
                    db_id = highest_id + 1
                # Save the design variables to the database
                req = "INSERT INTO solutions_" + SRM_name + " (id, "
                req += ", ".join(["angle_%i"%i for i in range(1,3)])
                # req += ", "
                # req += ", ".join(["TVC_y_%i"%i for i in range(1,6)])
                req += ", "
                req += ", ".join(["TVC_z_%i"%i for i in range(1,6)])
                req += ", "
                req += ", ".join(["spherical_motor_%i"%i for i in range(1,3)])
                req += ", "
                req += ", ".join(["%s_motor_%i"%(SRM_name, i) for i in range(1,len(SRM_1_geo_vars)+1)])
                req += ", dv_used) VALUES ("
                req += ", ".join(["?" for i in range(11+len(SRM_1_geo_vars))])
                req += ")"
                cur.execute(req, (db_id,)+tuple(dv)+(save_to_db,))

            # Set SRM thrust model from geometry models
            SRM_thrust_model_1 = SRM_thrust(SRM_1_model, A_t=0.065, epsilon=45)
            SRM_thrust_model_2 = SRM_thrust(SRM_2_model, A_t=0.005, epsilon=73, p_a=0)

            if plot_geo:
                ax = SRM_1_model.plot_geometry(add_title=False)
                ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
                ax.grid(False)
                ax.set(frame_on=False)
                plt.tight_layout()
                plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/motors/%i_1_geo.pdf"%(SRM_name, db_id))
                plt.close()
                ax = SRM_2_model.plot_geometry(add_title=False)
                ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
                ax.grid(False)
                ax.set(frame_on=False)
                plt.tight_layout()
                plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/motors/%i_2_geo.pdf"%(SRM_name, db_id))
                plt.close()
            if plot_thrust:
                SRM_thrust_model_1.simulate_full_burn(dt=2e-5, use_cpp=True)
                # M_P_s = [SRM_1_model.get_V_p()*SRM_thrust_model_1.rho_p]
                # for i in range(1,len(SRM_thrust_model_1.saved_burn_times)):
                #     M_P_s.append(M_P_s[-1] + SRM_thrust_model_1.m_dot_interpolator(SRM_thrust_model_1.saved_burn_times[i])*(SRM_thrust_model_1.saved_burn_times[i]-SRM_thrust_model_1.saved_burn_times[i-1]))
                # Get index of all magnitudes below 0
                ind_below_0 = [i for i in range(len(SRM_thrust_model_1.saved_magnitudes)) if SRM_thrust_model_1.saved_magnitudes[i] < 0]
                idx_0 = min(ind_below_0) if len(ind_below_0) != 0 else -1
                plt.figure(figsize=(4.5,3.25))
                plt.plot(SRM_thrust_model_1.saved_burn_times[:idx_0], np.asarray(SRM_thrust_model_1.saved_magnitudes[:idx_0])/1e3)
                plt.xlabel("Burn time [s]"), plt.ylabel("Thrust [kN]")
                plt.grid(), plt.tight_layout()
                plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/motors/%i_1_thrust.pdf"%(SRM_name, db_id))
                plt.close()
                print("Geo of %s plotted" % SRM_name)
                plt.figure(figsize=(4.5,3.25))
                SRM_thrust_model_2.simulate_full_burn(dt=1.5e-2)
                plt.plot(SRM_thrust_model_2.saved_burn_times, np.asarray(SRM_thrust_model_2.saved_magnitudes)/1e3)
                plt.xlabel("Burn time [s]"), plt.ylabel("Thrust [kN]")
                plt.grid(), plt.tight_layout()
                plt.savefig(sys.path[0]+"/plots/optimisation/results/%s/motors/%i_2_thrust.pdf"%(SRM_name, db_id))
                print("Thrust of %s plotted" % SRM_name)
                print("Burn time:", SRM_thrust_model_1.saved_burn_times[idx_0])
                mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
                mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p
                print("Launch mass of %.3f kg, second stage of %.3f kg" % (mass_1, mass_2))
                return 0
            
            if not skip_sim or better_accuracy:
                # inputs.append((launch_angle_1, launch_angle_2, TVC_angles_y, TVC_angles_z, SRM_thrust_model_1, SRM_thrust_model_2, False, db_id))
                inputs.append((launch_angle_1, launch_angle_2, 0, TVC_angles_z, SRM_thrust_model_1, SRM_thrust_model_2, False, db_id, better_accuracy))

        if save_to_db is not None:
            con.commit()
        con.close()

        # Get the fitness by running the simulations in parallel
        with MP.get_context("spawn").Pool(MP.cpu_count()//2) as pool:
            outputs = pool.starmap(MAV_ascent_sim, inputs)

        # Return the 1D list of fitnesses
        for i_res in range(n_dvs):
            if i_res in i_exist:
                h_p_score, h_a_score, mass_score = scores_exist.pop(0)
            else:
                h_p_score, h_a_score, mass_score = outputs.pop(0)
            fitnesses.append(h_p_score)
            fitnesses.append(h_a_score)
            fitnesses.append(mass_score)

        return fitnesses


if __name__ == "__main__":
    from thrust.models.multi_fin import multi_fin_SRM
    opti_problem = MAV_problem(0, multi_fin_SRM)