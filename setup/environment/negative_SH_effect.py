import enum
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
from matplotlib import pyplot as plt

# Tudatpy imports
from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import propagation_setup

from setup import ascent

models = ["PM", "SH D/O 2", "SH D/O 4", "SH D/O 8", "SH D/O 16"]
model_times = []

plt.figure(figsize=(10, 6))

for i, model in enumerate(models):
    MAV_ascent = ascent.MAV_ascent(
        launch_epoch = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17))),    # MAV-­LL­-01
        launch_lat = np.deg2rad(18.85),     # MAV­-LL-­03
        launch_lon = np.deg2rad(77.52),     # MAV­-LL-­03
        launch_h = -2.5e3,                  # MAV­-LL-­04
        mass_1 = [370, 185],                # MAV-­VM­-03 + LS p.10
        mass_2 = [80, 40],                  # LS p.10
        launch_angles = [np.deg2rad(45), np.deg2rad(90)],   # MAV­-LL-­06 + guesstimate # angle is w.r.t vertical
        thrust_models = [[9750, 293], [6750, 282]],         # magnitude, Isp for both stage, adapted from LS p.10
        target_orbit_h = 300e3,             # MAV­-OSO­-01
        target_orbit_i = np.deg2rad(25),    # MAV­-OSO­-03
        max_a = 15 * 9.80665,               # MAV­-LL-­02
        max_AoA = np.deg2rad(4)             # MAV­-LL-­05
    )

    MAV_ascent.create_bodies(stage=1, include_aero=False)
    accelerations_on_vehicle = MAV_ascent.create_accelerations(only_thrust_dict=True)
    if i == 0:
        acc = propagation_setup.acceleration.point_mass_gravity()
    else:
        acc = propagation_setup.acceleration.spherical_harmonic_gravity(2**i, 2**i)
    accelerations_on_vehicle["Mars"] = [acc]
    acceleration_dict = {MAV_ascent.current_name: accelerations_on_vehicle}
    MAV_ascent.acceleration_models = propagation_setup.create_acceleration_models(
        MAV_ascent.bodies, acceleration_dict, MAV_ascent.bodies_to_propagate, MAV_ascent.central_bodies
    )
    MAV_ascent.create_initial_state()
    MAV_ascent.create_dependent_variables_to_save(default=False)
    MAV_ascent.dependent_variables_to_save = [
        propagation_setup.dependent_variable.altitude( "MAV stage 1", "Mars" ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.thrust_acceleration_type, "MAV stage 1", "MAV stage 1")
    ]
    if i == 0:
        MAV_ascent.dependent_variables_to_save.append(
            propagation_setup.dependent_variable.single_acceleration_norm(
                propagation_setup.acceleration.point_mass_gravity_type, "MAV stage 1", "Mars")
        )
    else:
        MAV_ascent.dependent_variables_to_save.append(
            propagation_setup.dependent_variable.single_acceleration_norm(
                propagation_setup.acceleration.spherical_harmonic_gravity_type, "MAV stage 1", "Mars")
        )
    end_time = 80
    MAV_ascent.combined_termination_settings = propagation_setup.propagator.time_termination(MAV_ascent.initial_epoch + end_time)
    MAV_ascent.create_propagator_settings()
    MAV_ascent.create_integrator_settings(fixed_step=0.1)
    times, states, dep_vars = MAV_ascent.run_simulation()

    times = (times - times[0])/60
    altitudes = dep_vars[:,0]
    a_thrust = dep_vars[:,1]
    a_gravity = dep_vars[:,2]

    if i == 0:
        plt.plot(times, a_thrust, label="Thrust")
        idx_h_0 = np.argwhere(altitudes >= 0)[0][0]
        time_h_0 = times[idx_h_0]
        a_limits = [-2, 55]
        plt.vlines([time_h_0], a_limits[0], a_limits[1], color="k", linestyles="dotted")
        plt.ylim(a_limits[0], a_limits[1])
    plt.plot(times, a_gravity, label=model)

# plt.xlim(0, end_time/60)
plt.xlabel("Time [min]"), plt.ylabel("Acceleration [m/s$^2$]")
plt.legend()
plt.grid(), plt.tight_layout()

plt.show()