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
from matplotlib import pyplot as plt

# Tudatpy imports
# from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice

from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.anchor import anchor_SRM
from thrust.models.spherical import spherical_SRM
from thrust.models.tubular import tubular_SRM
# from thrust.solid_thrust import SRM_thrust
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust

use_SRM = True

# Max D=0.57m, most assumed in feasibility study was D=0.51m
if use_SRM:
    SRM_stage_1 = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.15)
    # SRM_stage_1 = rod_and_tube_SRM(R_o=0.25, R_mid=0.14, R_i=0.05, L=1.05)
    # SRM_stage_1 = anchor_SRM(R_o=0.26, R_i=0.165, N_a=6, w=0.025, r_f=0.005, delta_s=0.015, L=1.15)
    # SRM_stage_1 = tubular_SRM(R_o=0.25, R_i=0.16, L=1.15)
    SRM_thrust_model_1 = SRM_thrust(SRM_stage_1, A_t=0.065, epsilon=45)

    SRM_stage_2 = spherical_SRM(R_o=0.165, R_i=0.0875) # with multi-fin
    # SRM_stage_2 = spherical_SRM(R_o=0.19, R_i=0.0915) # with rod and tube
    # SRM_stage_2 = spherical_SRM(R_o=0.19, R_i=0.0915) # with anchor
    # SRM_stage_2 = spherical_SRM(R_o=0.19, R_i=0.12) # with tubular
    SRM_thrust_model_2 = SRM_thrust(SRM_stage_2, A_t=0.005, epsilon=73, p_a=0)
    # print("%.2f/207 kg of propellant"%SRM_thrust_model_1.M_p, "%.2f/29 kg of innert"%SRM_thrust_model_1.M_innert)
    # SRM_stage_1.plot_geometry()
    # plt.show()
    # print("%.2f/52 kg of propellant"%SRM_thrust_model_2.M_p, "%.2f/15 kg of innert"%SRM_thrust_model_2.M_innert)
    # SRM_stage_2.plot_geometry()
    # plt.show()

    mass_2 = 47.5 + SRM_thrust_model_2.M_innert + SRM_thrust_model_2.M_p
    mass_1 = 65 + mass_2 + SRM_thrust_model_1.M_innert + SRM_thrust_model_1.M_p
    print("Rocket mass of %.2f kg for section 1, %.2f kg for section 2." % (mass_1, mass_2))
else:
    mass_1, mass_2 = 370, 95
    # Magnitude, Isp, burn time
    SRM_thrust_model_1 = [11.5e3, 291, 54.5]    # Total impulse of 632.5e3 Ns
    SRM_thrust_model_2 = [4.5e3, 282, 22.5]     # Total impulse of 101.25e3 Ns

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
    launch_epoch = 0,#time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 2, 17))),    # MAV-­LL­-01
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
    body_fixed_thrust_direction_z=body_fixed_thrust_direction_z
)

# Setup and run simulation for both stages
stage_res = []
for stage in [1, 2]:
    MAV_ascent.create_bodies(stage=stage)
    MAV_ascent.create_accelerations(use_cpp=True)#thrust_fname=sys.path[0]+"/data/opti_thrust_%i.npz"%stage
    guidance_object = ascent_framework.FakeAeroGuidance()
    environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
    MAV_ascent.create_initial_state()
    MAV_ascent.create_dependent_variables_to_save()
    MAV_ascent.create_termination_settings(end_time=160*60)
    MAV_ascent.create_propagator_settings()
    MAV_ascent.create_integrator_settings(better_accuracy=True),#fixed_step=0.1)
    print("Running the simulation for stage %i" % stage)
    times, states, dep_vars = MAV_ascent.run_simulation()
    stage_res.append([times, states, dep_vars])
    final_h = max(dep_vars[:,1])
    print("Max altitude at end of stage %i at %.2f min of %.2f km, at %.2f km/s, vehicle mass of %.2f kg..." % \
        (stage, (times[-1]-MAV_ascent.launch_epoch)/60, final_h/1e3, dep_vars[-1,5]/1e3, states[-1,-1]))
    if stage == 1:
        t_b_1 = MAV_ascent.thrust.burn_time
        print("Stage 1 burn time = ", t_b_1)
        t_sep = times[-1]-MAV_ascent.launch_epoch
        idx_sep = len(times)-1
        if final_h < 0:
            break
    else:
        t_b_2 = MAV_ascent.thrust.burn_time
        print("Stage 2 burn time = ", t_b_2)
if stage == 1:
    # Extract results from first propagation only if stage 2 was not used
    times = stage_res[0][0]
    states = stage_res[0][1]
    dep_vars = stage_res[0][2]
else:
    # Combine results from both propagations
    times = np.concatenate((stage_res[0][0], stage_res[1][0]))
    states = np.concatenate((stage_res[0][1], stage_res[1][1]))
    dep_vars = np.concatenate((stage_res[0][2], stage_res[1][2]))

# np.savez(sys.path[0]+"/data/MAV_ascent_test.npz", times=times, states=states, dep_vars=dep_vars)

times = (times - times[0])/60

flight_path_angles = np.rad2deg(dep_vars[:,0])
altitudes = dep_vars[:,1]
force_coeffs = dep_vars[:,2:5]
airspeeds = dep_vars[:,5]
tot_accs = dep_vars[:,6]
mach_numbers = dep_vars[:,7]
mass = dep_vars[:,8]
angle_of_attacks = np.rad2deg(dep_vars[:,9])
positions = dep_vars[:,10:13]
a_SH = dep_vars[:,13]
a_thrust = dep_vars[:,14]
a_aero = dep_vars[:,15]
dyna_pressures = dep_vars[:,16]
velocities = dep_vars[:,17:20]
full_a_thrust = dep_vars[:,20:23]

print("Mach number ranges from %.2f to %.2f" % (min(mach_numbers), max(mach_numbers)))

final_cartesian_state = states[-1,:-1]

final_keplerian_state = element_conversion.cartesian_to_keplerian(
    final_cartesian_state,
    spice.get_body_gravitational_parameter("Mars")
)
final_a, final_e, final_i = final_keplerian_state[:3]
R_Mars = spice.get_average_radius("Mars")

h_peri, h_apo = final_a * (1 - final_e) - R_Mars, final_a * (1 + final_e) - R_Mars
print("Periapsis altitude = %.2f km, Apoapsis altitude = %.2f km" % (h_peri/1e3, h_apo/1e3))

exit()

try:
    idx_crop = np.where(times >= t_sep/60+3)[0][0]
except IndexError:
    idx_crop = -1

# idx_crop = np.where(times >= 20)[0][0]
dts = np.diff(times*60)
print("Minimum timestep was of %.3e s, maximum of %.3e s." % (np.ma.masked_array(dts, mask=dts==0).min(), max(dts)))

plt.figure(figsize=(9,5))
plt.plot(times[:idx_crop], dts[:idx_crop])
plt.grid()
plt.xlabel("Time [min]"), plt.ylabel("Time step [s]")
plt.yscale("log")
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/benchmark_dts.pdf")
plt.close()

# fig = plt.figure(figsize=(7, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.scatter(*positions[:idx_crop].T, color="C0", alpha=0.25)
# ax.scatter(*positions[0].T, color="C2")
# ax.scatter(*positions[idx_sep].T, color="C4")
# for i, position in enumerate(positions[:idx_crop]):
#     ax.plot(*np.array([position, position+velocities[i]*2.5e2]).T, color="C1", alpha=0.25)
#     ax.plot(*np.array([position, position+full_a_thrust[i]*5e3]).T, color="C3", alpha=0.25)
# ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
# plt.show()

import matplotlib
matplotlib.rcParams.update({'font.size': 16})

plt.figure(figsize=(7, 6.5))
plt.plot(times, altitudes/1e3)
plt.grid()
plt.xlabel("Time since launch [min]")
plt.ylabel("Altitude [km]")
plt.tight_layout()


plt.figure(figsize=(7, 6.5))
plt.plot(times[:idx_crop], mass[:idx_crop])
plt.grid()
plt.xlabel("Time since launch [min]")
plt.ylabel("Rocket mass [kg]")
plt.tight_layout()

plt.figure(figsize=(7, 6.5))
plt.plot(times[:idx_crop], tot_accs[:idx_crop], label="Total", linestyle="dotted", color="black")
plt.plot(times[:idx_crop], a_SH[:idx_crop], label="SH D/O 4")
plt.plot(times[:idx_crop], a_thrust[:idx_crop], label="Thrust")
plt.plot(times[:idx_crop], a_aero[:idx_crop], label="Aerodynamic")
plt.grid()
plt.xlabel("Time since launch [min]")
plt.ylabel("Acceleration [m/s$^2$]")
#ax5.set_yscale("log") # (uncomment this to see the distinction between the accelerations more clearly)
plt.legend()
plt.tight_layout()

plt.show()

# Create a figure with 5 subplots: a grid of 2x2, then one horizontal one at the bottom
fig = plt.figure(figsize=(14, 15))
gs = fig.add_gridspec(3, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 =fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax5 = fig.add_subplot(gs[2, :])

# Plot the altitude history
ax1.plot(times, altitudes/1e3)
ax1.grid()
ax1.set_xlabel("Time since launch [min]")
ax1.set_ylabel("Altitude [km]")

# Plot the airspeed history
ax2.plot(times, airspeeds)
ax2.grid()
ax2.set_xlabel("Time since launch [min]")
ax2.set_ylabel("Airspeed [m/s]")

# Plot the mass history in the first 10 minutes
ax3.plot(times[:idx_crop], mass[:idx_crop])
ax3.grid()
ax3.set_xlabel("Time since launch [min]")
ax3.set_ylabel("Rocket mass [kg]")

# Plot the dynamic pressure history in the first 10 minutes
ax4.plot(times[:idx_crop], dyna_pressures[:idx_crop]/1e3)
ymin, ymax = ax4.get_ylim()
xmin, xmax = ax4.get_xlim()
ax4.vlines(t_b_1/60, -1e3, ymax*1e3, color="k", linestyle="dotted")
ax4.hlines(3.5, -1e3, xmax*1e3, color="red", linestyle="dotted")
ax4.set_xlim(xmin, xmax)
ax4.set_ylim(ymin, ymax)
ax4.grid()
ax4.set_xlabel("Time since launch [min]")
ax4.set_ylabel("Dynamic pressure [kPa]")

# Plot the accelerations history in the first 10 minutes
ax5.plot(times[:idx_crop], tot_accs[:idx_crop], label="Total", linestyle="dotted", color="black")
ax5.plot(times[:idx_crop], a_SH[:idx_crop], label="SH D/O 4")
ax5.plot(times[:idx_crop], a_thrust[:idx_crop], label="Thrust")
ax5.plot(times[:idx_crop], a_aero[:idx_crop], label="Aerodynamic")
ax5.grid()
ax5.set_xlabel("Time since launch [min]")
ax5.set_ylabel("Acceleration [m/s$^2$]")
#ax5.set_yscale("log") # (uncomment this to see the distinction between the accelerations more clearly)
ax5.legend()

# Save some space using a tight layout, and show the figure
plt.tight_layout()
plt.show()


# plt.figure(figsize=(10, 6))
# plt.plot(times[:idx_crop], force_coeffs[:,0][:idx_crop], label="CD")
# plt.plot(times[:idx_crop], force_coeffs[:,2][:idx_crop], label="CL")
# plt.xlabel("Time [min]"), plt.ylabel("Force coefficient [-]")
# plt.legend()
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

# plt.show()