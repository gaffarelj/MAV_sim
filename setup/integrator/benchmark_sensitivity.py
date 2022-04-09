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
import glob
import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Tudatpy imports
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup, estimation_setup
from tudatpy import util

# Custom imports
from setup import ascent_framework_benchmarks
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust import SRM_thrust


# Main parameters
stage = 2                       # Stage for which to analyse the sensitivity analysis
powered = True                  # Is the stage be powered or not
allowable_errors = [15, 0.015]   # Allowable final error in position
deviation_exponents = np.linspace(-1, 0.35, 9, dtype=float)    # Exponent for the deviation from the nominal value


# Define initial mass for each ascent part
stage_initial_masses = [387.2554241260012, 187.22431444345327, 89.48046825001722, 61.428692147331795]

# Define factor by which to scale initial velocity error compared to position error (based on benchmark dt/error)
vel_factors = [10, 10, 100, 1000]

# If we are at lift-off, use the pre-defined initial state
if stage == 1 and powered:
    initial_state = None
# Load results from previous ascent part to get final state
else:
    if powered:
        f_name = "/data/best_integrator_dt/%i_X_dt_*.npz"%(stage-1)
        final_mass = stage_initial_masses[0] if stage == 1 else stage_initial_masses[2]
        vel_factor = vel_factors[0] if stage == 1 else vel_factors[2]
    else:
        f_name = "/data/best_integrator_dt/%i_V_dt_*.npz"%stage
        final_mass = stage_initial_masses[1] if stage == 1 else stage_initial_masses[3]
        vel_factor = vel_factors[1] if stage == 1 else vel_factors[3]
    previous_f_name = glob.glob(sys.path[0]+f_name)[0]
    previous_results = np.load(previous_f_name)
    # Get final (new initial) time, state, mass
    final_time, final_state = previous_results["times"][-1], previous_results["states"][-1][:6]
    initial_state = [final_time, final_state, final_mass]

# Get dt that was used for current ascent part from the saved filename
current_f_name = glob.glob(sys.path[0]+"/data/best_integrator_dt/%i_%s_dt_*.npz" % (stage, "V" if powered else "X"))[0]
current_dt = float(re.match(".+dt_([0-9]\.[0-9e+-]+).+", current_f_name).groups()[0])

# Define ascent
t0 = 0
body_fixed_thrust_direction_y = [
    [0, 0.05, 0.1, 0, 0.05],
    0
]
body_fixed_thrust_direction_z = [
    [0, -0.05, 0.0, 0.05, 0.05],
    0
]
SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
SRM_stage_2 = spherical_SRM(R_o=0.165, R_i=0.0915)
SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)
SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)
mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p
MAV_ascent = ascent_framework_benchmarks.MAV_ascent(
    launch_epoch = t0,
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
    body_fixed_thrust_direction_z=body_fixed_thrust_direction_z,
    powered=powered
)

# Setup simulation
MAV_ascent.dt = current_dt
MAV_ascent.create_bodies(stage=stage)
thrust_filename = glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0] if powered else None
MAV_ascent.create_accelerations(thrust_filename=thrust_filename)
guidance_object = ascent_framework_benchmarks.FakeAeroGuidance()
environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
MAV_ascent.create_initial_state(state=initial_state)
MAV_ascent.create_dependent_variables_to_save(default=False)
MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(MAV_ascent.current_name, "Mars"))
MAV_ascent.create_termination_settings(end_time=160*60)
MAV_ascent.create_integrator_settings(fixed_step=current_dt)

fig = plt.figure(figsize=(10, 7))
gs = fig.add_gridspec(2, 1)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
print("Propagating dynamics with dt = %.2e..." % current_dt)
if powered:
    deviation_exponents = np.insert(deviation_exponents, 0, -np.inf)
    baseline_states = None
    for i, dev in enumerate(deviation_exponents):
        initial_state_variation = [10**dev, 10**dev, 10**dev, 10**dev/vel_factor, 10**dev/vel_factor, 10**dev/vel_factor]
        print("Running with initial state deviation =", initial_state_variation)
        MAV_ascent.create_propagator_settings(initial_state_variation)
        states, dep_vars, f_evals = MAV_ascent.run_simulation(return_raw=True, return_count=True)
        # Save the baseline states
        if i == 0:
            baseline_states = util.result2array(states)
        else:
            states_array = util.result2array(states)
            states_diff = states_array - baseline_states
            # Plot the deviation in position
            ax1.plot(
                states_array[:,0]/60,
                np.linalg.norm(states_diff[:,1:4], axis=1),
                linestyle="dashed",
                label="%.2e m / %.2e m/s" % (10**dev, 10**dev/vel_factor)
            )
            # Plot the deviation in velocity
            ax2.plot(
                states_array[:,0]/60,
                np.linalg.norm(states_diff[:,4:], axis=1),
                linestyle="dashed",
                label="%.2e m / %.2e m/s" % (10**dev, 10**dev/vel_factor)
            )
else:
    MAV_ascent.create_propagator_settings()
    # Setup variational equations
    parameter_settings = estimation_setup.parameter.initial_states(MAV_ascent.propagator_settings, MAV_ascent.bodies)
    parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, MAV_ascent.bodies)

    # Propagate the dynamics
    with util.redirect_std():
        variational_equations_solver = numerical_simulation.SingleArcVariationalSimulator(
            MAV_ascent.bodies,
            MAV_ascent.integrator_settings,
            MAV_ascent.propagator_settings,
            parameters_to_estimate
        )

    # Extract simulation results
    states = variational_equations_solver.state_history
    dep_vars = variational_equations_solver.dynamics_simulator.dependent_variable_history
    state_transition_matrices = variational_equations_solver.state_transition_matrix_history

    # Compute state history with initial state deviation
    print("Computing sensitivity to deviations...")
    for dev in deviation_exponents:
        initial_state_variation = [10**dev, 10**dev, 10**dev, 10**dev/vel_factor, 10**dev/vel_factor, 10**dev/vel_factor]
        delta_initial_state_dict = dict()
        for epoch in state_transition_matrices:
            delta_initial_state_dict[epoch] = np.dot(state_transition_matrices[epoch], initial_state_variation)
        delta_initial_state_array = util.result2array(delta_initial_state_dict)

        # Plot the deviation in position
        ax1.plot(
            delta_initial_state_array[:,0]/60,
            np.linalg.norm(delta_initial_state_array[:,1:4], axis=1),
            linestyle="dashed",
            label="%.2e m / %.2e m/s" % (10**dev, 10**dev/vel_factor)
        )
        # Plot the deviation in velocity
        ax2.plot(
            delta_initial_state_array[:,0]/60,
            np.linalg.norm(delta_initial_state_array[:,4:], axis=1),
            linestyle="dashed",
            label="%.2e m / %.2e m/s" % (10**dev, 10**dev/vel_factor)
        )

# Plot allowable error
xlims = ax1.get_xlim()
ax1.hlines(allowable_errors[0], -1e3, 1e3, colors="orange")
ax1.set_xlim(xlims)
xlims = ax2.get_xlim()
ax2.hlines(allowable_errors[1], -1e3, 1e3, colors="orange")
ax2.set_xlim(xlims)

# Prettify
ax1.set_xlabel("Time [min]"), ax1.set_ylabel("Deviation in position [m]")
ax1.set_yscale("log")
ax1.grid()
# Add legend with allowable error
legend_elements = [Line2D([0], [0], color='orange', label="Allowable error (%.2e m)"%allowable_errors[0])]
lgd1 = ax1.legend(handles=legend_elements, loc='lower right')
# Add second legend above top subplot
ax1.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
           mode="expand", borderaxespad=0, ncol=3)
ax1.add_artist(lgd1)
ax2.set_xlabel("Time [min]"), ax2.set_ylabel("Deviation in velocity [m/s]")
ax2.set_yscale("log")
ax2.grid()
# Add legend with allowable error
legend_elements = [Line2D([0], [0], color='orange', label="Allowable error (%.2e m/s)"%allowable_errors[1])]
ax2.legend(handles=legend_elements, loc='lower right')
plt.tight_layout()
plt.savefig(sys.path[0] + "/plots/setup/integrator/sensitivity_%i_%s.pdf" % (stage, "V" if powered else "X"))
plt.show()