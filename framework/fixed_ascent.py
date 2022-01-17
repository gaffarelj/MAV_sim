import sys
# Add tudatpy path
sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Standard imports
import numpy as np
from datetime import datetime

# Tudatpy imports
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion, time_conversion
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment, environment_setup
from tudatpy.kernel.numerical_simulation import propagation, propagation_setup
from tudatpy.util import result2array

# Load the SPICE kernel
spice_interface.load_standard_kernels()

# Key parameters from requirements
launch_epoch = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17)))    # MAV-­LL­-01
launch_lat, launch_lon = np.deg2rad(18.85), np.deg2rad(77.52)   # MAV­-LL-­03
launch_h = -2.5e3                                               # MAV­-LL-­04
launch_mass = 400                                               # MAV-­VM­-03
launch_angle = np.deg2rad(90)                                   # MAV­-LL-­06
target_orbit_h = 300e3                                          # MAV­-OSO­-01
target_orbit_i = np.deg2rad(25)                                 # MAV­-OSO­-03
max_acceleration = 15 * 9.80665                                 # MAV­-LL-­02
max_AoA = np.deg2rad(4)                                         # MAV­-LL-­05

# Create bodies
bodies_to_create = ["Mars"]
body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "ECLIPJ2000")
body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")
bodies = environment_setup.create_system_of_bodies(body_settings)
central_bodies = ["Mars"]

# Create vehicle
bodies.create_empty_body("MAV")
bodies.get_body("MAV").set_constant_mass(launch_mass)
force_coefficients_files = {i: "data/coefficients/CF%s.dat" % l for i, l in enumerate(["x", "y", "z"])}
moment_coefficients_files = {i: "data/coefficients/CM%s.dat" % l for i, l in enumerate(["x", "y", "z"])}
coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_from_files(
    force_coefficient_files=force_coefficients_files,
    moment_coefficient_files=moment_coefficients_files,
    reference_length=0.43,
    reference_area=0.144,
    lateral_reference_length=0.43,
    moment_reference_point = [1.125, 0, 0],
    independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent]
)
environment_setup.add_aerodynamic_coefficient_interface(bodies, "MAV", coefficient_settings)
bodies_to_propagate = ["MAV"]

# Define accelerations
thrust_direction_settings = (
    propagation_setup.thrust.thrust_direction_from_state_guidance(
        central_body="Mars",
        is_colinear_with_velocity=True,
        direction_is_opposite_to_vector=False ) )
thrust_magnitude_settings = (
    propagation_setup.thrust.constant_thrust_magnitude(
        thrust_magnitude=9854, specific_impulse=293 ))
accelerations_on_rocket = dict(
    MAV=[
        propagation_setup.acceleration.thrust_from_direction_and_magnitude(
            thrust_direction_settings,
            thrust_magnitude_settings,
        )
    ],
    Mars=[
        # propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
        propagation_setup.acceleration.point_mass_gravity(),
        propagation_setup.acceleration.aerodynamic()
    ]
)

acceleration_dict = dict(MAV=accelerations_on_rocket)
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_dict, bodies_to_propagate, central_bodies
)

# Create initial state
initial_mars_fixed_state = element_conversion.spherical_to_cartesian_elementwise(
    radial_distance=bodies.get_body("Mars").shape_model.average_radius + launch_h,
    latitude=launch_lat,
    longitude=launch_lon,
    speed=0.1,
    flight_path_angle=launch_angle,
    heading_angle=0
)
initial_state = environment.transform_to_inertial_orientation(
    initial_mars_fixed_state, launch_epoch, bodies.get_body("Mars").rotation_model
)

# List dependent variables
dependent_variables_to_save = [
    propagation_setup.dependent_variable.flight_path_angle("MAV", "Mars"),
    propagation_setup.dependent_variable.altitude("MAV", "Mars"),
    propagation_setup.dependent_variable.aerodynamic_force_coefficients("MAV"),
    propagation_setup.dependent_variable.airspeed("MAV", "Mars"),
    propagation_setup.dependent_variable.total_acceleration_norm("MAV"),
    propagation_setup.dependent_variable.mach_number("MAV", "Mars"),
    propagation_setup.dependent_variable.body_mass("MAV"),
    propagation_setup.dependent_variable.angle_of_attack("MAV", "Mars"),
    propagation_setup.dependent_variable.relative_position("MAV", "Mars"),
    propagation_setup.dependent_variable.single_acceleration_norm(
        propagation_setup.acceleration.point_mass_gravity_type, "MAV", "Mars"),
    # propagation_setup.dependent_variable.single_acceleration_norm(
    #     propagation_setup.acceleration.spherical_harmonic_gravity_type, "MAV", "Mars"),
    propagation_setup.dependent_variable.single_acceleration_norm(
        propagation_setup.acceleration.thrust_acceleration_type, "MAV", "MAV"),
    propagation_setup.dependent_variable.single_acceleration_norm(
        propagation_setup.acceleration.aerodynamic_type , "MAV", "Mars")
]

# Create termination settings
termination_max_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
    dependent_variable_settings=propagation_setup.dependent_variable.altitude("MAV", "Mars"),
    limit_value=target_orbit_h,
    use_as_lower_limit=False)
termination_min_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
    dependent_variable_settings=propagation_setup.dependent_variable.altitude("MAV", "Mars"),
    limit_value=-5e3,
    use_as_lower_limit=True)
combined_termination_settings = propagation_setup.propagator.hybrid_termination(
[termination_max_altitude_settings, termination_min_altitude_settings], fulfill_single_condition=True)

# Create translational and mass propagators
translational_propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    combined_termination_settings,
    propagation_setup.propagator.encke,
    #propagation_setup.propagator.cowell,
    output_variables=dependent_variables_to_save
)
mass_rate_settings = dict(MAV=[propagation_setup.mass_rate.from_thrust()])
mass_rate_models = propagation_setup.create_mass_rate_models(
    bodies,
    mass_rate_settings,
    acceleration_models
)
mass_propagator_settings = propagation_setup.propagator.mass(
    bodies_to_propagate, mass_rate_models, [launch_mass], combined_termination_settings)
propagator_settings = propagation_setup.propagator.multitype(
    [translational_propagator_settings, mass_propagator_settings],
    combined_termination_settings,
    dependent_variables_to_save)

# Create integrator settings
initial_time_step = 1.0
minimum_time_step = 0.001
maximum_time_step = 500
tolerance = 1e-9
integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
    launch_epoch,
    initial_time_step,
    propagation_setup.integrator.rkf_78,
    minimum_time_step,
    maximum_time_step,
    relative_error_tolerance=tolerance,
    absolute_error_tolerance=tolerance)

# Propagate
dynamics_simulator = numerical_simulation.SingleArcSimulator(
    bodies, integrator_settings, propagator_settings, print_dependent_variable_data=True
)

# Analyse results
states = result2array(dynamics_simulator.state_history)
dep_vars = result2array(dynamics_simulator.dependent_variable_history)

times = dep_vars[:,0]
times = (times - times[0])/60
flight_path_angles = np.rad2deg(dep_vars[:,1])
altitudes = dep_vars[:,2]
force_coeffs = dep_vars[:,3:6]
airspeeds = dep_vars[:,6]
tot_accs = dep_vars[:,7]
mach_numbers = dep_vars[:,8]
mass = dep_vars[:,9]
angle_of_attacks = np.rad2deg(dep_vars[:,10])
a_pm = dep_vars[:,14]
a_thrust = dep_vars[:,15]
a_aero = dep_vars[:,16]

X, Y, Z = dep_vars[:,11], dep_vars[:,12], dep_vars[:,13]

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(times, a_pm, label="PM")
plt.plot(times, a_thrust, label="Thrust")
plt.plot(times, a_aero, label="Aero")
plt.xlabel("Time [min]"), plt.ylabel("Acceleration [m/s$^2$]")
plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, altitudes/1e3)
# plt.xlabel("Time [min]"), plt.ylabel("Altitude [km]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, airspeeds)
# plt.xlabel("Time [min]"), plt.ylabel("Airspeed [m/s]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, mass)
# plt.xlabel("Time [min]"), plt.ylabel("Mass [kg]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, force_coeffs[:,0])
# plt.xlabel("Time [min]"), plt.ylabel("CD [-]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, force_coeffs[:,2])
# plt.xlabel("Time [min]"), plt.ylabel("CL [-]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, tot_accs/9.80665)
# plt.xlabel("Time [min]"), plt.ylabel("g-load [-]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, flight_path_angles)
# plt.xlabel("Time [min]"), plt.ylabel("Flight path angle [deg]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(times, angle_of_attacks)
# plt.xlabel("Time [min]"), plt.ylabel("Angle of attack [deg]")
# plt.grid(), plt.tight_layout()

# plt.figure(figsize=(10, 6))
# plt.plot(np.sqrt(X**2 + Y**2)/1e3, Z/1e3)
# plt.xlabel("X-Y [km]"), plt.ylabel("Z [km]")
# plt.grid(), plt.tight_layout()

plt.show()