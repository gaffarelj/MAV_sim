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
from datetime import datetime
import copy
from scipy.interpolate import interp1d
import glob

# Tudatpy imports
from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.math import interpolators

# Custom imports
from setup import ascent_framework_segmented
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust import SRM_thrust

def resample(x, n=5000, kind='linear'):
    f = interp1d(np.linspace(0, 1, x.size), x, kind)
    return f(np.linspace(0, 1, n))

t0 = 0 #time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17)))
body_fixed_thrust_direction_y = [
    [0, 0.05, 0.1, 0, 0.05],
    0
]
body_fixed_thrust_direction_z = [
    [0, -0.05, 0.0, 0.05, 0.05],
    0
]
interpolator_settings = interpolators.lagrange_interpolation(8, boundary_interpolation=interpolators.use_boundary_value)

def run_all(dt, stage, powered=True, only_burn=False):

    SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
    SRM_stage_2 = spherical_SRM(R_o=0.165, R_i=0.0915)
    SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)
    SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)
    mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
    mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p

    # Define default ascent model
    MAV_ascent = ascent_framework_segmented.MAV_ascent(
        launch_epoch = t0,    # MAV-­LL­-01
        launch_lat = np.deg2rad(18.85),     # MAV­-LL-­03
        launch_lon = np.deg2rad(77.52),     # MAV­-LL-­03
        launch_h = -2.5e3,                  # MAV­-LL-­04
        mass_stages = [mass_1, mass_2],            
        launch_angles = [np.deg2rad(57.5), np.deg2rad(90)],       # MAV­-LL-­06 + guesstimate # angle is w.r.t vertical
        thrust_models = [SRM_thrust_model_1, SRM_thrust_model_2],
        target_orbit_h = 300e3,             # MAV­-OSO­-01
        target_orbit_i = np.deg2rad(25),    # MAV­-OSO­-03
        max_a = 15 * 9.80665,               # MAV­-LL-­02
        max_AoA = np.deg2rad(4),            # MAV­-LL-­05
        body_fixed_thrust_direction_y=body_fixed_thrust_direction_y,
        body_fixed_thrust_direction_z=body_fixed_thrust_direction_z,
        powered=powered
    )

    if only_burn:
        print("Runing stage %i burn sim with dt = %.3e" % (stage, dt))
        if stage == 1:
            times, magnitudes, *_, masses = SRM_thrust_model_1.simulate_full_burn(dt, make_interplators=False)
        elif stage == 2:
            times, magnitudes, *_, masses = SRM_thrust_model_2.simulate_full_burn(dt, make_interplators=False)

        fevals = len(times)
        times, magnitudes, masses = np.asarray(times), np.asarray(magnitudes), np.asarray(masses)

        times, magnitudes, masses = resample(times), resample(magnitudes), resample(masses)

        np.savez("setup/integrator/benchmark_sim_results/thrust_%i_dt_%.4e" % (stage, dt), \
            times=times, magnitudes=magnitudes, masses=masses)
        print("dt = %.3e s, stage = %i -> %.3e f evals" % (dt, stage, fevals))
    else:
        MAV_ascent.dt = dt
        # Setup and run simulation for stage 1 then 2
        print("Running with dt = %.3e s, stage = %i, %s" % (dt, stage, "powered" if powered else "unpowered"))
        MAV_ascent.create_bodies(stage=stage)
        thrust_filename = None
        if powered:
            thrust_filename = glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0]
        MAV_ascent.create_accelerations(thrust_filename=thrust_filename)
        guidance_object = ascent_framework_segmented.FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
        last_state_fname = None
        if not (powered and stage == 1):
            if not powered:
                f_name = "/data/best_integrator_dt/%i_V_dt_*.npz"%stage
            else:
                f_name = "/data/best_integrator_dt/%i_X_dt_*.npz"%(stage-1)
            last_state_fname = glob.glob(sys.path[0]+f_name)[0]
        MAV_ascent.create_initial_state(last_state_fname)
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(MAV_ascent.current_name, "Mars"))
        MAV_ascent.create_termination_settings(end_time=160*60)
        MAV_ascent.create_propagator_settings()
        MAV_ascent.create_integrator_settings(fixed_step=dt)
        states, dep_vars, f_evals = MAV_ascent.run_simulation(return_raw=True, return_count=True)
        times = np.asarray(list(states.keys()))
        if len(times) <= 8:
            print("Not enough time steps made with dt = %.4e (only %i)..." % (dt, len(times)))
            np.savez("setup/integrator/benchmark_sim_results/%i_%s_dt_%.4e" % (stage, "V" if powered else "X", dt), times=None, states=None, dep_vars=None, f_evals=None)
        else:
            # Resample times
            times = resample(times)
            # Interpolate state and dep vars based on resampled times
            states_interpolator = interpolators.create_one_dimensional_vector_interpolator(states, interpolator_settings)
            dep_vars_interpolator = interpolators.create_one_dimensional_vector_interpolator(dep_vars, interpolator_settings)
            states = np.asarray([states_interpolator.interpolate(epoch) for epoch in times])
            dep_vars = np.asarray([dep_vars_interpolator.interpolate(epoch) for epoch in times])

            np.savez("setup/integrator/benchmark_sim_results/%i_%s_dt_%.4e" % (stage, "V" if powered else "X", dt), times=times, states=states, dep_vars=dep_vars, f_evals=f_evals)
            print("dt = %.3e s, stage = %i, %s -> %.3e f evals" % (dt, stage, "powered" if powered else "unpowered", f_evals))