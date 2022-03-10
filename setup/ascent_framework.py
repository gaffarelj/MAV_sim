import numpy as np
from scipy import interpolate
import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment, environment_setup
from tudatpy.kernel.numerical_simulation import propagation, propagation_setup
from tudatpy.util import result2array
from tudatpy.kernel.math import root_finders

from thrust.MAV_thrust import MAV_thrust

# Load the SPICE kernel
spice_interface.load_standard_kernels()

class FakeAeroGuidance(propagation.AerodynamicGuidance):

    def __init__(self):
        # Call the base class constructor
        propagation.AerodynamicGuidance.__init__(self)

    def updateGuidance(self, current_time):
        # Update angle of attack
        self.angle_of_attack = np.deg2rad(1.5) * np.sin(current_time*np.pi/750)

class MAV_ascent:

    def __init__(self, launch_epoch, launch_lat, launch_lon, launch_h, mass_stages, \
        launch_angles, thrust_models, target_orbit_h, target_orbit_i, max_a, max_AoA, body_fixed_thrust_direction):
        self.launch_epoch = launch_epoch
        self.launch_lat, self.launch_lon = launch_lat, launch_lon
        self.launch_h = launch_h
        self.stage_1_wet_mass = mass_stages[0]
        self.stage_2_wet_mass = mass_stages[1]
        self.launch_angles = launch_angles
        self.thrust_models = thrust_models
        self.target_orbit_h = target_orbit_h
        self.target_orbit_i = target_orbit_i
        self.max_acceleration = max_a
        self.max_angle_of_attack = max_AoA
        self.last_h = -np.inf
        self.body_fixed_thrust_direction = body_fixed_thrust_direction

    def create_bodies(self, stage, include_aero=True):
        self.current_stage = stage
        self.current_name = "MAV stage %i" % self.current_stage
        # Create Mars
        bodies_to_create = ["Mars"]
        # TODO: investigate why changing reference frame orientation seems to change results (thrust orientation?)
        body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "IAU_Mars")
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
            CDs = CDs = [1.5986, 1.7730, 1.7659, 1.7397, 1.7166, 1.7022, 1.6930, 1.6871, 1.6843, 1.6830, 1.6837, 1.6843, 1.6815]
            hs =   np.array([100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500])*1e3
            cs = interpolate.interp1d(hs, CDs, kind="quadratic", bounds_error=False, fill_value="extrapolate")
            def get_CD(dep_vars):
                h = dep_vars[0]
                return [cs(h), 0, 0]
            coefficient_settings = environment_setup.aerodynamic_coefficients.custom(
                get_CD,
                reference_area=0.144,
                independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent]
            )
        if include_aero:
            environment_setup.add_aerodynamic_coefficient_interface(self.bodies, self.current_name, coefficient_settings)
        self.bodies_to_propagate = [self.current_name]

    def create_accelerations(self, only_thrust_dict=False):
        # Define thrust
        self.thrust = MAV_thrust(
            self,
            self.launch_angles[self.current_stage-1],
            self.thrust_models[self.current_stage-1],
            self.body_fixed_thrust_direction[self.current_stage-1])
        thrust_direction_settings = propagation_setup.thrust.custom_thrust_direction(self.thrust.get_thrust_orientation)
        thrust_magnitude_settings = propagation_setup.thrust.custom_thrust_magnitude(
            self.thrust.get_thrust_magnitude,
            self.thrust.get_specific_impulse,
            self.thrust.is_thrust_on)

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
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
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

    def create_dependent_variables_to_save(self, default=True):
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
        else:
            self.dependent_variables_to_save = []

    def create_termination_settings(self, end_time=200*60):

        # def is_vehicle_falling(_time):
        #     dh = self.current_body.flight_conditions.altitude - self.last_h
        #     self.last_h = self.current_body.flight_conditions.altitude
        #     position = self.current_body.position
        #     velocity = self.current_body.velocity
        #     unit_vector_pos = position / np.linalg.norm(position)
        #     unit_vector_vel = velocity / np.linalg.norm(velocity)
        #     dot_product = np.dot(unit_vector_vel, unit_vector_pos)
        #     angle = np.pi/2-np.arccos(dot_product)
        #     # print(dh, np.rad2deg(angle))
        #     # if dh < 0:
        #     #     input()
        #     return np.rad2deg(angle) < 2# and self.last_h > 1e3
        #     return dh < 0

        if self.current_stage == 1:
            termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                limit_value=-10e3,
                use_as_lower_limit=True)
            # termination_apogee_settings = propagation_setup.propagator.custom_termination(is_vehicle_falling)
            termination_apogee_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.flight_path_angle(self.current_name, "Mars"),
                limit_value=0,
                use_as_lower_limit=True,
                terminate_exactly_on_final_condition=True,
                termination_root_finder_settings=root_finders.secant(
                    maximum_iteration=5,
                    maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result)
                )
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(
            [termination_apogee_settings, termination_min_altitude_settings], fulfill_single_condition=True)
        else:
            termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings=propagation_setup.dependent_variable.altitude(self.current_name, "Mars"),
                limit_value=50e3,
                use_as_lower_limit=True)
            termination_max_time_settings = propagation_setup.propagator.time_termination(self.initial_epoch + end_time)
            self.combined_termination_settings = propagation_setup.propagator.hybrid_termination(
            [termination_max_time_settings, termination_min_altitude_settings], fulfill_single_condition=True)


    def create_propagator_settings(self):
        translational_propagator_settings = propagation_setup.propagator.translational(
            self.central_bodies,
            self.acceleration_models,
            self.bodies_to_propagate,
            self.initial_state,
            self.combined_termination_settings,
            propagation_setup.propagator.cowell,
            output_variables=self.dependent_variables_to_save
        )
        
        if self.thrust.thrust_type == "constant":
            mass_rate_settings = {self.current_name:[propagation_setup.mass_rate.from_thrust()]}
        else:
            mass_rate_settings = {self.current_name:[propagation_setup.mass_rate.custom(self.thrust.get_mass_flow)]}
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
            coefficients = propagation_setup.integrator.rkf_78
        else:
            initial_time_step = 1e-2
            minimum_time_step = 1e-8
            maximum_time_step = 500
            tolerance = 1e-18
            coefficients = propagation_setup.integrator.rkf_78
        self.integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
            self.initial_epoch,
            initial_time_step,
            coefficients,
            minimum_time_step,
            maximum_time_step,
            relative_error_tolerance=tolerance,
            absolute_error_tolerance=tolerance)

    def run_simulation(self, return_raw=False, return_count=False):
        dynamics_simulator = numerical_simulation.SingleArcSimulator(
            self.bodies,
            self.integrator_settings,
            self.propagator_settings,
            print_dependent_variable_data=False,
            print_state_data=False
        )

        raw_states, raw_dep_vars = dynamics_simulator.state_history, dynamics_simulator.dependent_variable_history
        states = result2array(raw_states)
        dep_vars = result2array(raw_dep_vars)
        self.times = dep_vars[:,0]
        self.states = states[:,1:]
        self.dep_vars = dep_vars[:,1:]
        ret  = None
        if return_raw:
            ret = raw_states, raw_dep_vars
        else:
            ret = self.times, self.states, self.dep_vars
        if return_count:
            f_evals = list(dynamics_simulator.cumulative_number_of_function_evaluations.values())[-1]
            return tuple((*ret, f_evals))
        return ret