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
from tudatpy.kernel.astro import time_conversion

import ascent


MAV_ascent = ascent.MAV_ascent(
    launch_epoch = 0,#time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17)))    # MAV-­LL­-01
    launch_lat = np.deg2rad(18.85),     # MAV­-LL-­03
    launch_lon = np.deg2rad(77.52),     # MAV­-LL-­03
    launch_h = -2.5e3,                  # MAV­-LL-­04
    mass_1 = [400, 185],                # MAV-­VM­-03 + LS p.10
    mass_2 = [100, 45],                 # LS p.10
    launch_angle = np.deg2rad(45),      # MAV­-LL-­06
    target_orbit_h = 300e3,             # MAV­-OSO­-01
    target_orbit_i = np.deg2rad(25),    # MAV­-OSO­-03
    max_a = 15 * 9.80665,               # MAV­-LL-­02
    max_AoA = np.deg2rad(4),            # MAV­-LL-­05
    staging_altitude = 200e3
)

# Setup and run simulation for stage 1
stage_res = []
for stage in [1, 2]:
    print("Running the simulation for stage %i" % stage)
    MAV_ascent.create_bodies(stage=stage)
    MAV_ascent.create_accelerations()
    MAV_ascent.create_initial_state()
    MAV_ascent.create_dependent_variables_to_save()
    MAV_ascent.create_termination_settings()
    MAV_ascent.create_propagator_settings()
    MAV_ascent.create_integrator_settings()
    times, states, dep_vars = MAV_ascent.run_simulation()
    stage_res.append([times, states, dep_vars])
    if stage==1:
        print("Altitude at second stage start:", dep_vars[-1,1])

# Combine results from both propagations
times = np.concatenate((stage_res[0][0], stage_res[1][0]))
states = np.concatenate((stage_res[0][1], stage_res[1][1]))
dep_vars = np.concatenate((stage_res[0][2], stage_res[1][2]))

times = (times - times[0])/60

flight_path_angles = np.rad2deg(dep_vars[:,0])
altitudes = dep_vars[:,1]
force_coeffs = dep_vars[:,2:5]
airspeeds = dep_vars[:,5]
tot_accs = dep_vars[:,6]
mach_numbers = dep_vars[:,7]
mass = dep_vars[:,8]
angle_of_attacks = np.rad2deg(dep_vars[:,9])
a_pm = dep_vars[:,13]
a_thrust = dep_vars[:,14]
a_aero = dep_vars[:,15]


X, Y, Z = dep_vars[:,10], dep_vars[:,11], dep_vars[:,12]

import matplotlib.pyplot as plt

idx_crop = np.where(times >= 12)[0][0]

plt.figure(figsize=(10, 6))
plt.plot(times[:idx_crop], a_pm[:idx_crop], label="SH")
plt.plot(times[:idx_crop], a_thrust[:idx_crop], label="Thrust")
plt.plot(times[:idx_crop], a_aero[:idx_crop], label="Aero")
plt.plot(times[:idx_crop], tot_accs[:idx_crop], label="Total", linestyle="dotted")
plt.xlabel("Time [min]"), plt.ylabel("Acceleration [m/s$^2$]")
plt.legend()
plt.grid(), plt.tight_layout()

plt.figure(figsize=(10, 6))
plt.plot(times, altitudes/1e3)
plt.xlabel("Time [min]"), plt.ylabel("Altitude [km]")
plt.grid(), plt.tight_layout()

plt.figure(figsize=(10, 6))
plt.plot(times, airspeeds)
plt.xlabel("Time [min]"), plt.ylabel("Airspeed [m/s]")
plt.grid(), plt.tight_layout()

plt.figure(figsize=(10, 6))
plt.plot(times[:idx_crop], mass[:idx_crop])
plt.xlabel("Time [min]"), plt.ylabel("Mass [kg]")
plt.grid(), plt.tight_layout()

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