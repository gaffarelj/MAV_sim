import numpy as np
import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment, environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.util import result2array
from tudatpy.kernel.math import root_finders

from thrust import MAV_thrust

# Load the SPICE kernel
spice_interface.load_standard_kernels()

class MAV_ascent:

    def __init__(self, launch_epoch, launch_lat, launch_lon, launch_h, mass_1, mass_2, \
        launch_angles, thrust_magnitudes, target_orbit_h, target_orbit_i, max_a, max_AoA, staging_altitude):
        self.launch_epoch = launch_epoch
        self.launch_lat, self.launch_lon = launch_lat, launch_lon
        self.launch_h = launch_h
        self.stage_1_wet_mass, self.stage_1_dry_mass = mass_1[0], mass_1[1]
        self.stage_2_wet_mass, self.stage_2_dry_mass = mass_2[0], mass_2[1]
        self.launch_angles = launch_angles
        self.thrust_magnitudes = thrust_magnitudes
        self.target_orbit_h = target_orbit_h
        self.target_orbit_i = target_orbit_i
        self.max_acceleration = max_a
        self.max_angle_of_attack = max_AoA
        self.staging_altitude = staging_altitude

    def create_bodies(self, stage):
        self.current_stage = stage
        self.current_name = "MAV stage %i" % self.current_stage
        # Create Mars
        bodies_to_create = ["Mars"]
        body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "ECLIPJ2000")
        body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        self.central_bodies = ["Mars"]

        # Create vehicle
        self.bodies.create_empty_body(self.current_name)
        self.current_body = self.bodies.get(self.current_name)
        if stage == 1:
            self.current_body.set_constant_mass(self.stage_1_wet_mass)
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
            coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
                0.144,
                constant_force_coefficient=[1.5, 0, 0.3]
            )
        environment_setup.add_aerodynamic_coefficient_interface(self.bodies, self.current_name, coefficient_settings)
        self.bodies_to_propagate = [self.current_name]

    def create_accelerations(self, only_thrust_dict=False):
        # Define thrust
        if self.current_stage == 1:
            thrust = MAV_thrust(self, self.launch_angles[0], self.thrust_magnitudes[0], 293)
        else:
            thrust = MAV_thrust(self, self.launch_angles[1], self.thrust_magnitudes[1], 282)
        #thrust_direction_settings = propagation_setup.thrust.custom_thrust_orientation(thrust.get_thrust_orientation)
        thrust_direction_settings = propagation_setup.thrust.custom_thrust_direction(thrust.get_thrust_orientation)
        thrust_magnitude_settings = propagation_setup.thrust.custom_thrust_magnitude(thrust.get_thrust_magnitude, thrust.get_specific_impulse, thrust.is_thrust_on)

        if only_thrust_dict:
            accelerations_on_vehicle = {
                self.current_name: [
                    propagation_setup.acceleration.thrust_from_direction_and_magnitude(
                        thrust_direction_settings,
                        thrust_magnitude_settings
                    )
                ]
            }
            return accelerations_on_vehicle

        # Define accelerations
        accelerations_on_vehicle = {
            self.current_name: [
                propagation_setup.acceleration.thrust_from_direction_and_magnitude(
                    thrust_direction_settings,
                    thrust_magnitude_settings
                )
            ],
            "Mars": [
                propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
                propagation_setup.acceleration.aerodynamic()
            ]
        }

        acceleration_dict = {self.current_name: accelerations_on_vehicle}
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies, acceleration_dict, self.bodies_to_propagate, self.central_bodies
        )

    def create_initial_state(self):
        if self.current_stage == 1:
            self.initial_epoch = self.launch_epoch
            initial_mars_fixed_state = element_conversion.spherical_to_cartesian_elementwise(
                radial_distance=self.bodies.get_body("Mars").shape_model.average_radius + self.launch_h,
                latitude=self.launch_lat,
                longitude=self.launch_lon,
                speed=1.0,
                flight_path_angle=self.launch_angles[0],
                heading_angle=0
            )
            self.initial_state = environment.transform_to_inertial_orientation(
                initial_mars_fixed_state, self.initial_epoch, self.bodies.get_body("Mars").rotation_model
            )
        elif self.current_stage == 2:
            self.initial_epoch = self.times[-1]
            try:
                last_state = self.states[-1]
            except AttributeError:
                raise ValueError("The last states of the first stage could not be found. Please make sure that the simulation has been run for the first stage before attempting to run for the second stage.")
            self.initial_state = last_state[:6]

    def create_dependent_variables_to_save(self, only_a=None):
        if only_a is not None:
            self.dependent_variables_to_save = []
        else:
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
                    propagation_setup.acceleration.aerodynamic_type, self.current_name, "Mars")
            ]

    def create_termination_settings(self):
        termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
            limit_value=-10e3,
            use_as_lower_limit=True)
        if self.current_stage == 1:
            termination_max_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                limit_value=self.staging_altitude,
                use_as_lower_limit=False,
                terminate_exactly_on_final_condition=True,
                termination_root_finder_settings=root_finders.secant(
                    maximum_iteration=3,
                    maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result))
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(
            [termination_max_altitude_settings, termination_min_altitude_settings], fulfill_single_condition=True)
        else:
            termination_max_time_settings = propagation_setup.propagator.time_termination(self.initial_epoch + 200*60)
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(
            [termination_max_time_settings, termination_min_altitude_settings], fulfill_single_condition=True)


    def create_propagator_settings(self):
        translational_propagator_settings = propagation_setup.propagator.translational(
            self.central_bodies,
            self.acceleration_models,
            self.bodies_to_propagate,
            self.initial_state,
            self.combined_termination_settings,
            propagation_setup.propagator.gauss_keplerian,
            output_variables=self.dependent_variables_to_save
        )
        
        mass_rate_settings = {self.current_name:[propagation_setup.mass_rate.from_thrust()]}
        mass_rate_models = propagation_setup.create_mass_rate_models(
            self.bodies,
            mass_rate_settings,
            self.acceleration_models
        )
        if self.current_stage == 1:
            launch_mass = self.stage_1_wet_mass
        else:
            launch_mass = self.stage_2_wet_mass
        mass_propagator_settings = propagation_setup.propagator.mass(
            self.bodies_to_propagate, mass_rate_models, [launch_mass], self.combined_termination_settings)
        self.propagator_settings = propagation_setup.propagator.multitype(
            [translational_propagator_settings, mass_propagator_settings],
            self.combined_termination_settings,
            self.dependent_variables_to_save)

    def create_integrator_settings(self, fixed_step=None):
        if fixed_step is not None:
            initial_time_step = fixed_step
            minimum_time_step = fixed_step
            maximum_time_step = fixed_step
            tolerance = 1
            coefficients = propagation_setup.integrator.rkf_45
        else:
            initial_time_step = 1.0
            minimum_time_step = 0.0001
            maximum_time_step = 500
            tolerance = 1e-14
            coefficients = propagation_setup.integrator.rkf_78
        self.integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
            self.initial_epoch,
            initial_time_step,
            coefficients,
            minimum_time_step,
            maximum_time_step,
            relative_error_tolerance=tolerance,
            absolute_error_tolerance=tolerance)

    def run_simulation(self, return_raw=False):
        dynamics_simulator = numerical_simulation.SingleArcSimulator(
            self.bodies, self.integrator_settings, self.propagator_settings, print_dependent_variable_data=False
        )

        raw_states, raw_dep_vars = dynamics_simulator.state_history, dynamics_simulator.dependent_variable_history
        states = result2array(raw_states)
        dep_vars = result2array(raw_dep_vars)
        self.times = dep_vars[:,0]
        self.states = states[:,1:]
        self.dep_vars = dep_vars[:,1:]
        if return_raw:
            return raw_states, raw_dep_vars
        return self.times, self.states, self.dep_vars