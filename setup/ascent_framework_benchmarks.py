# Add tudatpy path
import sys, os
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
from scipy import interpolate

# Tudatpy imports
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice
from tudatpy.kernel.numerical_simulation import environment, environment_setup
from tudatpy.kernel.numerical_simulation import propagation, propagation_setup
from tudatpy.util import result2array
from tudatpy.kernel.math import root_finders

# Custom imports
from thrust.MAV_thrust import MAV_thrust

class FakeAeroGuidance(propagation.AerodynamicGuidance):
    """
    Class for aerodynamic guidance, varying the angle of attack over time from -1.5deg to 1.5deg.
    """

    def __init__(self):
        # Call the base class constructor
        propagation.AerodynamicGuidance.__init__(self)

    def updateGuidance(self, current_time):
        # Update angle of attack
        self.angle_of_attack = np.deg2rad(1.5) * np.sin(current_time*np.pi/750)

class MAV_ascent:
    """
    Mars Ascent Vehicle ascent class. Used to setup the ascent simulation of the two stages.
    """
    def __init__(self, launch_epoch, launch_lat, launch_lon, launch_h, mass_stages, \
        launch_angles, thrust_models, target_orbit_h, target_orbit_i, max_a, max_AoA, \
        body_fixed_thrust_direction_y, body_fixed_thrust_direction_z, powered):

        # Save launch values
        self.launch_epoch = launch_epoch            # Launch epoch in seconds since J2000
        self.launch_lon = launch_lon                # Launch longitude in radians
        self.launch_lat = launch_lat                # Launch latitude in radians
        self.launch_h = launch_h                    # Launch altitude in meters
        self.stage_1_wet_mass = mass_stages[0]      # First section wet mass (stage 1 + 2) in kg
        self.stage_2_wet_mass = mass_stages[1]      # Second stage wet mass in kg
        self.launch_angles = launch_angles          # Angle at which first and second stage ignition starts, in radians

        # Save the thrust models used for both stages
        self.thrust_models = thrust_models          # Either a list with [magnitude [N], Isp [s], burn time [s]] or a thrust.solid_thrust.SRM_thrust (for both stages)

        # Save the target orbit altitude and inclination
        self.target_orbit_h = target_orbit_h        # Target altitude in meters
        self.target_orbit_i = target_orbit_i        # Target inclination in radians

        # Save the ascent constraints
        self.max_acceleration = max_a               # Maximum acceleration in m/s2
        self.max_angle_of_attack = max_AoA          # Maximum angle of attack in radians

        # Save the nozzle deflection in both planes (TVC)
        self.body_fixed_thrust_direction_y = body_fixed_thrust_direction_y  # List of floats or float representing the deflections in radians, one for each stage
        self.body_fixed_thrust_direction_z = body_fixed_thrust_direction_z  # Same

        self.powered = powered
        self.dt = 0.01

        # Load the SPICE kernel
        spice.load_standard_kernels()

    def create_bodies(self, stage, include_aero=True):
        # Save which stage is now running, and create its name
        self.current_stage = stage
        self.current_name = "MAV stage %i" % self.current_stage

        # Create Mars body with exponential atmosphere
        bodies_to_create = ["Mars"]
        body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "J2000")
        body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        self.central_bodies = ["Mars"]

        # Create the rocket stage vehicle
        self.bodies.create_empty_body(self.current_name)
        self.current_body = self.bodies.get(self.current_name)

        # Apply different aerodynamic coefficients and initial mass for both stages
        if stage == 1:
            self.current_body.set_constant_mass(self.stage_1_wet_mass)
            # Use the force and moment aerodynamic coefficients from Missile DATCOM
            force_coefficients_files = {i: sys.path[0] + "/data/coefficients/CF%s.dat" % l for i, l in enumerate(["x", "y", "z"])}
            moment_coefficients_files = {i: sys.path[0] + "/data/coefficients/CM%s.dat" % l for i, l in enumerate(["x", "y", "z"])}
            coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_from_files(
                force_coefficient_files=force_coefficients_files,
                moment_coefficient_files=moment_coefficients_files,
                reference_length=0.43,
                reference_area=0.144,
                lateral_reference_length=0.43,
                moment_reference_point = [1.125, 0, 0],
                independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent]
            )
        else:
            self.current_body.set_constant_mass(self.stage_2_wet_mass)
            # Use the drag coefficients from DSMC
            CDs = CDs = [1.5986, 1.7730, 1.7659, 1.7397, 1.7166, 1.7022, 1.6930, 1.6871, 1.6843, 1.6830, 1.6837, 1.6843, 1.6815]
            hs =   np.array([100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500])*1e3
            cs = interpolate.interp1d(hs, CDs, kind="quadratic", bounds_error=False, fill_value="extrapolate")
            def get_CD(dep_vars):
                # Function that returns the linearly interpolated drag coefficient as a function of altitude
                h = dep_vars[0]
                return [cs(h), 0, 0]
            coefficient_settings = environment_setup.aerodynamic_coefficients.custom(
                get_CD,
                reference_area=0.144,
                independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent]
            )
        # Only include the aerodynamic coefficients if specified (otherwise the mass will still be set)
        if include_aero:
            environment_setup.add_aerodynamic_coefficient_interface(self.bodies, self.current_name, coefficient_settings)

        # Set the rocket stage as the body to propagate
        self.bodies_to_propagate = [self.current_name]

    def create_accelerations(self, only_thrust_dict=False, thrust_filename=None):
        if self.powered:
            # Setup the MAV thrust class from the thrust models input to this ascent class
            self.thrust = MAV_thrust(
                self,
                self.launch_angles[self.current_stage-1],
                self.thrust_models[self.current_stage-1],
                self.body_fixed_thrust_direction_y[self.current_stage-1],
                self.body_fixed_thrust_direction_z[self.current_stage-1],
                dt=self.dt,
                thrust_filename=thrust_filename)

            # Define the thrust acceleration direction and magnitude from the thrust class
            thrust_direction_settings = propagation_setup.thrust.custom_thrust_direction(self.thrust.get_thrust_orientation)
            thrust_magnitude_settings = propagation_setup.thrust.custom_thrust_magnitude(
                self.thrust.get_thrust_magnitude,
                self.thrust.get_specific_impulse,
                self.thrust.is_thrust_on)

            # Define the thrust acceleration of the vehicle on itself using the direction and magnitude
            accelerations_on_vehicle = {
                    self.current_name: [
                        propagation_setup.acceleration.thrust_from_direction_and_magnitude(
                            thrust_direction_settings,
                            thrust_magnitude_settings
                        )
                    ]
                }
        else:
            accelerations_on_vehicle = {}

        # Add environmental accelerations (Mars SH up to D/O 4 and aerodynamics)
        accelerations_on_vehicle["Mars"] = [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic()
            ]

        # Compile the dict of accelerations acting on the vehicle, and create the models for them
        acceleration_dict = {self.current_name: accelerations_on_vehicle}
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies, acceleration_dict, self.bodies_to_propagate, self.central_bodies
        )

    def create_initial_state(self, filename=None, state=None):
        if filename is not None:
            saved_results = np.load(filename)
            final_state = saved_results["states"][-1]
            self.initial_epoch = saved_results["times"][-1]
            # Set the body initial mass if the vehicle is un-powered (otherwise keep wet mass that was set in create_bodies)
            if not self.powered: 
                self.current_body.set_constant_mass(final_state[-1])
            # Set the second stage initial state as the last one from the first stage (ignoring mass)
            self.initial_state = final_state[:6]
        elif state is not None:
            self.initial_epoch = state[0]
            self.initial_state = state[1]
            self.current_body.set_constant_mass(state[2])
        # If we are with the first stage, use the defined initial state
        elif self.current_stage == 1:
            # Set the initial epoch as the user-specified launch epoch
            self.initial_epoch = self.launch_epoch

            # Define the initial state in the Mars fixed reference frame, from spherical elements
            initial_mars_fixed_state = element_conversion.spherical_to_cartesian_elementwise(
                radial_distance=self.bodies.get_body("Mars").shape_model.average_radius + self.launch_h,
                latitude=self.launch_lat,
                longitude=self.launch_lon,
                speed=1.0,
                flight_path_angle=self.launch_angles[0],
                heading_angle=0
            )

            # Convert the initial state to the inertial frame
            self.initial_state = environment.transform_to_inertial_orientation(
                initial_mars_fixed_state, self.initial_epoch, self.bodies.get_body("Mars").rotation_model
            )
        # If we are with the second stage, re-use the final state of the first stage
        elif self.current_stage == 2:
            # Set the initial epoch as the last one from the first stage propagation
            self.initial_epoch = self.times[-1]

            # Try to access the last state from the first stage propagation
            try:
                last_state = self.states[-1]
            except AttributeError:
                raise ValueError("The last states of the first stage could not be found. Please make sure that the simulation has been run for the first stage before attempting to run for the second stage.")

            # Set the second stage initial state as the last one from the first stage (ignoring mass)
            self.initial_state = last_state[:6]

    def create_dependent_variables_to_save(self, default=True):
        # Use a pre-defined set of dependent variables (their names should be self-explanatory)
        self.dependent_variables_to_save = []

    def create_termination_settings(self, end_time=200*60):
        # Terminations:
        # * Stage 1:
        #   - powered: burnoff (time = burn_time)
        #   - unpowered: apogee
        # * Stage 2:
        #   - powered: burnoff
        #   - unpowered: end_time

        # For the first stage, terminate at apogee or below a certain altitude
        if self.powered:
            self.combined_termination_settings = propagation_setup.propagator.time_termination(
                self.initial_epoch + self.thrust.burn_time,
                terminate_exactly_on_final_condition=True)

        elif self.current_stage == 1:
            self.combined_termination_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.flight_path_angle(self.current_name, "Mars"),
            limit_value=0,
            use_as_lower_limit=True,
            terminate_exactly_on_final_condition=True,
            termination_root_finder_settings=root_finders.secant(
                maximum_iteration=5,
                maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result)
            )
        elif self.current_stage == 2:
            self.combined_termination_settings = propagation_setup.propagator.time_termination(self.launch_epoch + end_time)

    def create_propagator_settings(self):
        # Create translational propagator settings, using the defined bodies/models
        translational_propagator_settings = propagation_setup.propagator.translational(
            self.central_bodies,
            self.acceleration_models,
            self.bodies_to_propagate,
            self.initial_state,
            self.combined_termination_settings,
            propagation_setup.propagator.cowell,
            output_variables=self.dependent_variables_to_save
        )

        if self.powered:
        
            # If the thrust model is setup for a constant thrust, define mass rate settings in relation to the thrust
            if self.thrust.thrust_type == "constant":
                mass_rate_settings = {self.current_name:[propagation_setup.mass_rate.from_thrust()]}
            # Otherwise, use the mass flow computed from the thrust model
            else:
                mass_rate_settings = {self.current_name:[propagation_setup.mass_rate.custom(self.thrust.get_mass_flow)]}
            
            # Create the mass rate model from the mass rate settings
            mass_rate_models = propagation_setup.create_mass_rate_models(
                self.bodies,
                mass_rate_settings,
                self.acceleration_models
            )

            # Set the initial mass
            if self.current_stage == 1:
                launch_mass = self.stage_1_wet_mass
            else:
                launch_mass = self.stage_2_wet_mass

            # Define the mass propagator settings
            mass_propagator_settings = propagation_setup.propagator.mass(
                self.bodies_to_propagate, mass_rate_models, [launch_mass], self.combined_termination_settings)

            # Define combined propagator settings including the translation and mass propagators
            self.propagator_settings = propagation_setup.propagator.multitype(
                [translational_propagator_settings, mass_propagator_settings],
                self.combined_termination_settings,
                self.dependent_variables_to_save)
        else:
            self.propagator_settings = translational_propagator_settings

    def create_integrator_settings(self, fixed_step=None):
        # If specified, setup a fixed step integrator by setting the tolerance to infinity
        if fixed_step is not None:
            initial_time_step = fixed_step
            minimum_time_step = fixed_step
            maximum_time_step = fixed_step
            tolerance = np.inf
        # Otherwise, define the variable step integrator with a tolerance of 1e-18
        else:
            initial_time_step = 1e-3
            minimum_time_step = 1e-8
            maximum_time_step = 60
            tolerance = 1e-18
        # Use RKF7(8) coefficients
        coefficients = propagation_setup.integrator.rkf_78
        self.integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
            self.initial_epoch,
            initial_time_step,
            coefficients,
            minimum_time_step,
            maximum_time_step,
            relative_error_tolerance=tolerance,
            absolute_error_tolerance=tolerance,
            maximum_factor_increase=2,
            minimum_factor_increase=0.05)

    def run_simulation(self, return_raw=False, return_count=False):
        # Run the ascent simulation (do not print the state or dependent variable content)
        dynamics_simulator = numerical_simulation.SingleArcSimulator(
            self.bodies,
            self.integrator_settings,
            self.propagator_settings,
            print_dependent_variable_data=False,
            print_state_data=False
        )

        # Get the dictionaries with the state and dependent variable histories
        raw_states, raw_dep_vars = dynamics_simulator.state_history, dynamics_simulator.dependent_variable_history
        
        # Convert the results to numpy arrays (always convert since the times will be re-used)
        states = result2array(raw_states)
        dep_vars = result2array(raw_dep_vars)
        self.times = dep_vars[:,0]
        self.states = states[:,1:]
        self.dep_vars = dep_vars[:,1:]
        
        # Return the results...
        ret  = None
        # ... as dictionaries
        if return_raw:
            ret = raw_states, raw_dep_vars
        # ... as numpy arrays
        else:
            ret = self.times, self.states, self.dep_vars
        # ... including the number of function evaluations
        if return_count:
            f_evals = list(dynamics_simulator.cumulative_number_of_function_evaluations.values())[-1]
            return tuple((*ret, f_evals))
        return ret