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
import sqlite3

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy import util
from tudatpy.kernel.math import interpolators

# Custom imports
from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust import SRM_thrust


def run_ascent(method, coefficients, *args):

    # Define fixed and variable step integrator
    def fixed_step_integrator_settings(init_time, coefficients, dt):
        return propagation_setup.integrator.runge_kutta_fixed_step_size(
                init_time, dt, coefficients
        )

    def variable_step_integrator_settings(init_time, coefficients, tolerance):
        initial_time_step = 1e-4
        minimum_time_step = 1e-8
        maximum_time_step = 60
        return propagation_setup.integrator.runge_kutta_variable_step_size(
            init_time,
            initial_time_step,
            coefficients,
            minimum_time_step,
            maximum_time_step,
            relative_error_tolerance=tolerance,
            absolute_error_tolerance=tolerance,
            maximum_factor_increase=2,
            minimum_factor_increase=0.05)

    SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
    SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)

    SRM_stage_2 = spherical_SRM(R_o=0.165, R_i=0.0915)
    SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)

    mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
    mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p

    # Max deflection of 5 deg
    body_fixed_thrust_direction_y = [
        [0, 0.05, 0.1, 0, 0.05],
        0   # second stage has no TVC
    ]
    body_fixed_thrust_direction_z = [
        [0, -0.05, 0.0, 0.05, 0.05],
        0   # second stage has no TVC
    ]

    MAV_ascent = ascent_framework.MAV_ascent(
        launch_epoch = 0,
        launch_lat = np.deg2rad(18.85),
        launch_lon = np.deg2rad(77.52),
        launch_h = -2.5e3,
        mass_stages = [mass_1, mass_2],
        launch_angles = [np.deg2rad(57.5), np.deg2rad(90)],
        thrust_models = [SRM_thrust_model_1, SRM_thrust_model_2],
        target_orbit_h = 300e3,
        target_orbit_i = np.deg2rad(25),
        max_a = 15 * 9.80665,
        max_AoA = np.deg2rad(4),
        body_fixed_thrust_direction_y=body_fixed_thrust_direction_y,
        body_fixed_thrust_direction_z=body_fixed_thrust_direction_z
    )

    # Setup and run simulation for both stages
    times, states, f_evals, t_sep = [], [], 0, 0
    print("Running %s %s %s" % (method, str(coefficients).split(".")[-1], args))
    for stage in [1, 2]:
        MAV_ascent.create_bodies(stage=stage)
        MAV_ascent.create_accelerations()
        guidance_object = ascent_framework.FakeAeroGuidance()
        environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
        MAV_ascent.create_initial_state()
        MAV_ascent.create_dependent_variables_to_save(default=False)
        MAV_ascent.create_termination_settings(end_time=160*60)
        MAV_ascent.create_propagator_settings()
        MAV_ascent.create_integrator_settings(empty_settings=True)
        # Set integrator settings
        if method == "fixed":
            dt = args[0]
            integrator_settings = fixed_step_integrator_settings(MAV_ascent.initial_epoch, coefficients, dt)
        else:
            tolerance = args[0]
            integrator_settings = variable_step_integrator_settings(MAV_ascent.initial_epoch, coefficients, tolerance)
        MAV_ascent.integrator_settings = integrator_settings
        times_i, states_i, _dep_vars_i, f_evals_i = MAV_ascent.run_simulation(return_count=True)
        f_evals += f_evals_i
        times.extend(times_i), states.extend(states_i)
        if stage == 1:
            t_sep = times[-1]

    # if times[-1] < 140*60:
    #     print("Simulation ended after %.2f minutes." % (times[-1]/60))
    #     return None

    # Convert the simulation results to a dict
    sim_dict = {epoch: state for epoch, state in zip(times, states)}

    # Load the benchmark results
    benchmark = np.load(sys.path[0]+"/data/best_integrator_dt/full_benchmark.npz")
    times_benchmark = benchmark["times"]
    states_benchmark = benchmark["states"]

    # Convert the benchmark results to a dict
    benchmark_dict = {epoch: state for epoch, state in zip(times_benchmark, states_benchmark)}

    # Compute the difference between the simulation and the benchmark
    interpolator_settings = interpolators.linear_interpolation()
    baseline_results_interpolator = interpolators.create_one_dimensional_vector_interpolator(benchmark_dict, interpolator_settings)
    new_results_interpolator = interpolators.create_one_dimensional_vector_interpolator(sim_dict, interpolator_settings)
    diff_times_stage_1 = np.linspace(0, t_sep-35, 500) # stop 35 before the separation because benchmark data is truncated
    states_diff_stage_1 = np.asarray([new_results_interpolator.interpolate(epoch) - baseline_results_interpolator.interpolate(epoch) for epoch in diff_times_stage_1])
    diff_times_stage_2 = np.linspace(t_sep+25, times[-1]-60, 500)
    states_diff_stage_2 = np.asarray([new_results_interpolator.interpolate(epoch) - baseline_results_interpolator.interpolate(epoch) for epoch in diff_times_stage_2])
    diff_times = np.concatenate((diff_times_stage_1, diff_times_stage_2))
    states_diff = np.concatenate((states_diff_stage_1, states_diff_stage_2))

    # Compute maximum error in position and velocity
    pos_error = np.fabs(np.linalg.norm(states_diff[:,0:3], axis=1))
    final_pos_error = pos_error[-1]
    vel_error = np.fabs(np.linalg.norm(states_diff[:,3:6], axis=1))
    final_vel_error = vel_error[-1]
    print("For %s %s %s, final error at %.2f min in position of %.2e [m] / in velocity of %.2e [m/s], with %.2e f evals (vs %.2e)" % \
        (method, str(coefficients).split(".")[-1], args, times[-1]/60, final_pos_error, final_vel_error, f_evals, benchmark["f_evals"]))

    # Connect to the database
    con = sqlite3.connect(sys.path[0]+"/setup/integrator/integrator_tuning/database.db")
    cur = con.cursor()

    res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=? AND dt=?", ("fixed", str(coefficients).split(".")[-1], dt))
    if res.fetchone() is None:
        cur.execute("INSERT INTO integrator (method, coefficients, end_time, f_evals, error_pos, error_vel, dt) VALUES (?, ?, ?, ?, ?, ?, ?)", ("fixed", str(coefficients).split(".")[-1], times[-1], f_evals, final_pos_error, final_vel_error, dt))
    else:
        cur.execute("UPDATE integrator SET end_time=?, f_evals=?, error_pos=?, error_vel=? WHERE method=? AND coefficients=? AND dt=?", (times[-1], f_evals, final_pos_error, final_vel_error, "fixed", str(coefficients).split(".")[-1], dt))
    
    con.commit()
    con.close()

    if False:
        from matplotlib import pyplot as plt
        # Plot state error over time
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(diff_times/60, pos_error, label="Position [m]")
        ax.plot(diff_times/60, vel_error, label="Velocity [m/s]")
        ax.plot(diff_times/60, np.fabs(states_diff[:,6]), label="Mass [kg]")
        ax.set_xlabel("Time [min]"), ax.set_ylabel("Absolute state error norm")
        ax.set_yscale("log")
        ax.grid()
        ax.legend()
        plt.tight_layout()

        # Plot time step used over time
        dts = np.diff(times)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(np.asarray(times[:-1])/60, dts)
        ax.set_xlabel("Time [min]"), ax.set_ylabel("Time step [s]")
        ax.grid()
        ax.set_yscale("log")
        plt.tight_layout()
        plt.show()