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
        body_fixed_thrust_direction_y, body_fixed_thrust_direction_z):

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

        # Load the SPICE kernel
        spice.load_standard_kernels()

    def create_bodies(
        self,
        stage,
        include_aero=True,
        include_radiation_pressure=False,
        more_bodies=False,
        custom_exponential_model=False,
        use_MCD=False,
        use_GRAM=False,
        use_new_coeffs=False,
        add_sun=False
    ):
        # Save which stage is now running, and create its name
        self.current_stage = stage
        self.current_name = "MAV stage %i" % self.current_stage

        # Create Mars body with exponential atmosphere
        bodies_to_create = ["Mars"]
        if include_radiation_pressure or add_sun:
            bodies_to_create = ["Mars", "Sun"]
        elif more_bodies:
            bodies_to_create = ["Mars", "Sun", "Jupiter", "Venus", "Earth", "Saturn"]
        body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "J2000")
        if custom_exponential_model:
            def rho_expo(h):
                if h <= 25e3:
                    rho_0 = 0.0159
                    h_s = 11049
                else:
                    rho_0 = 0.0525
                    h_s = 7295
                return rho_0 * np.exp(-h/h_s)
            body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature(
                rho_expo,
                constant_temperature=215,
                specific_gas_constant=197,
                ratio_of_specific_heats=1.3
            )
        elif use_MCD or use_GRAM:
            # Use the GRAM (2010) atmospheric model
            if use_GRAM:
                from atmosphere import GRAM
                GRAM_rho = GRAM.GRAM_atmo()
                get_density = GRAM_rho.get_density
            elif use_MCD:
                from atmosphere import MCD
                mcd = MCD.mcd_interface()
                get_density = mcd.density
            body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_four_dimensional_constant_temperature(
                get_density,
                constant_temperature=215,
                specific_gas_constant=197,
                ratio_of_specific_heats=1.3
            )
        else:
            body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        self.central_bodies = ["Mars"]

        # Create the rocket stage vehicle
        self.bodies.create_empty_body(self.current_name)
        self.current_body = self.bodies.get(self.current_name)

        # Apply different aerodynamic coefficients and initial mass for both stages
        if stage == 1:
            self.current_body.set_constant_mass(self.stage_1_wet_mass)
            add_txt = ""
            if use_new_coeffs:
                add_txt = "new_"
            # Use the force and moment aerodynamic coefficients from Missile DATCOM
            force_coefficients_files = {i: sys.path[0] + "/data/coefficients/%sCF%s.dat" % (add_txt, l) for i, l in enumerate(["x", "y", "z"])}
            moment_coefficients_files = {i: sys.path[0] + "/data/coefficients/%sCM%s.dat" % (add_txt, l) for i, l in enumerate(["x", "y", "z"])}
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

        if include_radiation_pressure:
            # Add solar radiation settings
            reference_area_radiation = 2.8*0.57 if self.current_stage == 1 else 0.57*2.8*0.25
            radiation_pressure_coefficient = 1.2
            occulting_bodies = ["Mars"]
            radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
                "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
            environment_setup.add_radiation_pressure_interface(self.bodies, self.current_name, radiation_pressure_settings)

        # Set the rocket stage as the body to propagate
        self.bodies_to_propagate = [self.current_name]

    def create_accelerations(
        self,
        only_thrust_dict=False,
        thrust_fname=None,
        thrust_dt=None,
        use_cpp=False,
        better_precision=False
    ):
        # Setup the MAV thrust class from the thrust models input to this ascent class
        self.thrust = MAV_thrust(
            self,
            self.launch_angles[self.current_stage-1],
            self.thrust_models[self.current_stage-1],
            self.body_fixed_thrust_direction_y[self.current_stage-1],
            self.body_fixed_thrust_direction_z[self.current_stage-1],
            thrust_filename=thrust_fname,
            dt=thrust_dt,
            use_cpp=use_cpp)

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

        # Return the acceleration if only thrust is to be included
        if only_thrust_dict:
            return accelerations_on_vehicle

        # Add environmental accelerations (Mars SH up to D/O 4 and aerodynamics)
        SH_order = 6 if better_precision else 4
        accelerations_on_vehicle["Mars"] = [
                propagation_setup.acceleration.spherical_harmonic_gravity(SH_order, SH_order),
                propagation_setup.acceleration.aerodynamic()
            ]
        if better_precision:
            accelerations_on_vehicle["Sun"] = [
                propagation_setup.acceleration.point_mass_gravity()
            ]

        # Compile the dict of accelerations acting on the vehicle, and create the models for them
        acceleration_dict = {self.current_name: accelerations_on_vehicle}
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies, acceleration_dict, self.bodies_to_propagate, self.central_bodies
        )

    def create_initial_state(self):
        # If we are with the first stage, use the defined initial state
        if self.current_stage == 1:
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
        if default:
            self.dependent_variables_to_save = [
                propagation_setup.dependent_variable.flight_path_angle(self.current_name, "Mars"),
                propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                propagation_setup.dependent_variable.aerodynamic_force_coefficients(self.current_name),
                propagation_setup.dependent_variable.airspeed(self.current_name, "Mars"),
                propagation_setup.dependent_variable.total_acceleration_norm(self.current_name),
                propagation_setup.dependent_variable.mach_number(self.current_name, "Mars"),
                propagation_setup.dependent_variable.body_mass(self.current_name),
                propagation_setup.dependent_variable.angle_of_attack(self.current_name, "Mars"),
                propagation_setup.dependent_variable.relative_position(self.current_name, "Mars"),
                propagation_setup.dependent_variable.single_acceleration_norm(
                    propagation_setup.acceleration.spherical_harmonic_gravity_type, self.current_name, "Mars"),
                propagation_setup.dependent_variable.single_acceleration_norm(
                    propagation_setup.acceleration.thrust_acceleration_type, self.current_name, self.current_name),
                propagation_setup.dependent_variable.single_acceleration_norm(
                    propagation_setup.acceleration.aerodynamic_type, self.current_name, "Mars"),
                propagation_setup.dependent_variable.dynamic_pressure( self.current_name ),
                propagation_setup.dependent_variable.relative_velocity(self.current_name, "Mars"),
                propagation_setup.dependent_variable.single_acceleration(
                    propagation_setup.acceleration.thrust_acceleration_type, self.current_name, self.current_name)
            ]
        # Otherwise, set the dependent variables to save as an empty list
        else:
            self.dependent_variables_to_save = []

    def create_termination_settings(self, end_time=200*60, exact_time=False, print_progress=False, cpu_time_termination=None):

        self.i = 0
        def print_t(t):
            if self.i % 500 == 0:
                altitude = self.current_body.flight_conditions.altitude
                print("Time %.2f s: altitude = %.2f km" % (t, altitude/1e3), end="\r")
            self.i += 1
            return False
        fake_term_settings = propagation_setup.propagator.custom_termination(print_t)

        if cpu_time_termination is not None:
            cpu_time_termination_settings = propagation_setup.propagator.cpu_time_termination( cpu_time_termination )

        termination_abs_max_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
            limit_value=5000e3,
            use_as_lower_limit=False)

        # For the first stage, terminate at apogee or below a certain altitude
        if self.current_stage == 1:
            # Define termination settings to finish when the vehicle get below -10km
            termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                limit_value=-10e3,
                use_as_lower_limit=True)
            
            # Define termination settings to finish exactly at apogee, when the flight path becomes 0
            termination_apogee_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.flight_path_angle(self.current_name, "Mars"),
                limit_value=0,
                use_as_lower_limit=True,
                terminate_exactly_on_final_condition=True,
                termination_root_finder_settings=root_finders.secant(
                    maximum_iteration=5,
                    maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result)
                )
            
            # Combine both termination settings, to finish as soon as one of the two is reached
            list_terminations_settings = [termination_min_altitude_settings, termination_apogee_settings, termination_abs_max_altitude_settings]
            if cpu_time_termination is not None:
                list_terminations_settings.append( cpu_time_termination_settings )
            if print_progress:
                list_terminations_settings.append( fake_term_settings )
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(list_terminations_settings, fulfill_single_condition=True)
        # For the second stage, terminate after a given time or below a given altitude
        else:
            # Terminate if the altitude becomes lower than 25km
            termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                limit_value=25e3,
                use_as_lower_limit=True)

            # Terminate after a given time has elapsed since the launch epoch
            termination_max_time_settings = propagation_setup.propagator.time_termination(
                self.launch_epoch + end_time,
                terminate_exactly_on_final_condition=exact_time)
            
            # Combine both termination settings, to finish as soon as one of the two is reached
            list_terminations_settings = [termination_max_time_settings, termination_min_altitude_settings, termination_abs_max_altitude_settings]
            if cpu_time_termination is not None:
                list_terminations_settings.append( cpu_time_termination_settings )
            if print_progress:
                list_terminations_settings.append( fake_term_settings )
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(list_terminations_settings, fulfill_single_condition=True)

    def create_propagator_settings(self, propagator_type=propagation_setup.propagator.cowell):
        # Create translational propagator settings, using the defined bodies/models
        translational_propagator_settings = propagation_setup.propagator.translational(
            self.central_bodies,
            self.acceleration_models,
            self.bodies_to_propagate,
            self.initial_state,
            self.combined_termination_settings,
            propagator_type,
            output_variables=self.dependent_variables_to_save
        )
        
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

    def create_integrator_settings(self, fixed_step=None, empty_settings=False, better_accuracy=False):
        if empty_settings:
            self.integrator_settings = None
        else:
            # If specified, setup a fixed step integrator by setting the tolerance to infinity
            if fixed_step is not None:
                initial_time_step = fixed_step
                minimum_time_step = fixed_step
                maximum_time_step = fixed_step
                tolerance = np.inf
            # Otherwise, define the variable step integrator with a tolerance of 1e-18
            else:
                initial_time_step = 4.5e-5
                minimum_time_step = 3e-5
                maximum_time_step = 60
                tolerance = 1e-7
                if better_accuracy:
                    tolerance = 1e-13
            # Use RKF4(5) coefficients
            coefficients = propagation_setup.integrator.rkf_45
            self.integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
                self.initial_epoch,
                initial_time_step,
                coefficients,
                minimum_time_step,
                maximum_time_step,
                relative_error_tolerance=tolerance,
                absolute_error_tolerance=tolerance,
                maximum_factor_increase=1.01,#+0.1,
                minimum_factor_increase=0.01,#+0.1,
                throw_exception_if_minimum_step_exceeded=False,
                assess_termination_on_minor_steps=True)

    def run_simulation(self, return_raw=False, return_count=False, return_success_status=False):
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
            if return_success_status:
                return tuple((*ret, f_evals, dynamics_simulator.integration_completed_successfully))
            return tuple((*ret, f_evals))
        return ret