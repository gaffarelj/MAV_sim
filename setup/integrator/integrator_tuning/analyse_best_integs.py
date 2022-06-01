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
from matplotlib import pyplot as plt
import time as T
import glob

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import element_conversion

# Custom imports
from setup import ascent_framework
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust

# Parameters
plot_thrust = False

# Define problem
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

timings = []

# Setup and run simulation for both stages
stage_res = []
for stage in [1, 2]:
    MAV_ascent.create_bodies(stage=stage)
    t0 = T.time()
    MAV_ascent.create_accelerations(use_cpp=(stage==1))#, thrust_fname=glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0])
    t1 = T.time()
    guidance_object = ascent_framework.FakeAeroGuidance()
    environment_setup.set_aerodynamic_guidance(guidance_object, MAV_ascent.current_body, silence_warnings=True)
    MAV_ascent.create_initial_state()
    MAV_ascent.create_dependent_variables_to_save(default=False)
    MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.altitude(MAV_ascent.current_name, "Mars"))
    MAV_ascent.dependent_variables_to_save.append(propagation_setup.dependent_variable.body_mass(MAV_ascent.current_name))
    MAV_ascent.create_termination_settings(end_time=160*60, exact_time=True)
    MAV_ascent.create_propagator_settings()
    MAV_ascent.create_integrator_settings()
    print("Running the simulation for stage %i" % stage)
    times, states, dep_vars = MAV_ascent.run_simulation()
    t2 = T.time()
    timings.append(t1-t0), timings.append(t2-t1)
    stage_res.append([times, states, dep_vars])
    final_h = max(dep_vars[:,0])
    print("Max altitude at end of stage %i at %.2f min of %.2f km, vehicle mass of %.2f kg..." % \
        (stage, (times[-1]-MAV_ascent.launch_epoch)/60, final_h/1e3, states[-1,-1]))
    if stage == 1:
        t_b_1 = MAV_ascent.thrust.burn_time
        t_sep = times[-1]-MAV_ascent.launch_epoch
        idx_sep = len(times)-1
        if final_h < 0:
            break
    else:
        t_b_2 = MAV_ascent.thrust.burn_time

    if plot_thrust:
        print("Plotting thrust...")
        # Get thrust used during simulation
        thrust_times = np.arange(0, MAV_ascent.thrust.burn_time, 0.1)
        thrust_magnitudes = [MAV_ascent.thrust.get_thrust_magnitude(t, use_time_input=True)/1e3 for t in thrust_times]
        thrust_mass_flows = [MAV_ascent.thrust.get_mass_flow(t, use_time_input=True) for t in thrust_times]
        fig1 = plt.figure(figsize=(9, 5))
        plt.plot(thrust_times, thrust_magnitudes, label="Simulated")
        fig2 = plt.figure(figsize=(9, 5))
        plt.plot(thrust_times, thrust_mass_flows, label="Simulated")
        # Get benchmark thrust
        thrust_fname = glob.glob(sys.path[0]+"/data/best_integrator_dt/thrust_%i_dt_*.npz"%stage)[0]
        MAV_ascent.thrust.thrust_model.simulate_full_burn(filename=thrust_fname)
        MAV_ascent.thrust.magnitude_function = MAV_ascent.thrust.thrust_model.magnitude_interpolator
        MAV_ascent.thrust.m_dot_function = MAV_ascent.thrust.thrust_model.m_dot_interpolator
        thrust_times_benchmark = np.arange(0, MAV_ascent.thrust.thrust_model.saved_burn_times[-1], 0.1)
        thrust_magnitudes_benchmark = [MAV_ascent.thrust.get_thrust_magnitude(t, use_time_input=True)/1e3 for t in thrust_times_benchmark]
        thrust_mass_flows_benchmark = [MAV_ascent.thrust.get_mass_flow(t, use_time_input=True) for t in thrust_times_benchmark]
        plt.figure(fig1)
        plt.plot(thrust_times_benchmark, thrust_magnitudes_benchmark, linestyle="dotted", label="Benchmark")
        plt.xlabel("Time [s]")
        plt.ylabel("Thrust [kN]")
        plt.title("Thrust profile for stage %i" % stage)
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.figure(fig2)
        plt.plot(thrust_times_benchmark, thrust_mass_flows_benchmark, linestyle="dotted", label="Benchmark")
        plt.xlabel("Time [s]")
        plt.ylabel("Mass flow [kg/s]")
        plt.grid()
        plt.legend()
        plt.tight_layout()

        plt.show()

print("CPU time repartition:", timings)
    
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

# Plot results
fig = plt.figure(figsize=(9, 5))
gs = fig.add_gridspec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
# Only plot mass up to 13min
try:
    idx_crop = np.where(times >= 20*60)[0][0]
except IndexError:
    idx_crop = -1
ax1.plot(times/60, dep_vars[:,0]/1e3, label="Tuned integrator")
ax2.plot(times[:idx_crop]/60, dep_vars[:idx_crop,-1], label="Tuned integrator")

# Load and plot benchmark results
f_name = "/data/best_integrator_dt/full_benchmark.npz"
saved_results = np.load(sys.path[0]+f_name)
bench_times = saved_results["times"]
print(bench_times[-1]/60, times[-1]/60)
bench_states = saved_results["states"]
bench_dep_vars = saved_results["dep_vars"]
bench_altitudes = bench_dep_vars[:, 0]
bench_masses = bench_states[:, -1]
idx_crop_benchmark = np.where(bench_times >= 20*60)[0][0]
ax1.plot(bench_times/60, bench_altitudes/1e3, linestyle="dotted", label="Benchmark")
ax2.plot(bench_times[:idx_crop_benchmark]/60, bench_masses[:idx_crop_benchmark], linestyle="dotted", label="Benchmark")

# Compute difference in final state
final_state_diff = states[-1] - bench_states[-1]
pos_diff = np.linalg.norm(final_state_diff[:3])
vel_diff = np.linalg.norm(final_state_diff[3:6])
print("Position difference: %.2f km" % (pos_diff/1e3))
print("Velocity difference: %.2f m/s" % (vel_diff))

# Prettify plot
ax1.set_xlabel("Time [min]")
ax1.set_ylabel("Altitude [km]")
ax2.set_xlabel("Time [min]")
ax2.set_ylabel("Mass [kg]")
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.show()

# Plot dt over time
fig = plt.figure(figsize=(9, 5))
dts = list(np.diff(times))
dts.append(dts[-1])
plt.plot(times[:idx_crop]/60, dts[:idx_crop])
plt.xlabel("Time [min]")
plt.ylabel("Time step [s]")
plt.yscale("log")
plt.grid()
plt.tight_layout()
plt.show()

# Plot apoapsis and periapsis over time
mu_Mars = MAV_ascent.bodies.get_body("Mars").gravitational_parameter
r_Mars = MAV_ascent.bodies.get_body("Mars").shape_model.average_radius
inclinaisons, h_a_s, h_p_s = [], [], []
idx_min = len(states[idx_sep:])//6
for state in states[idx_sep+idx_min:]:
    kep_elements = element_conversion.cartesian_to_keplerian(state[:6], mu_Mars)
    a, e = kep_elements[:2]
    i = np.rad2deg(kep_elements[2])
    h_a = a*(1-e)-r_Mars
    h_p = a*(1+e)-r_Mars
    inclinaisons.append(i)
    h_a_s.append(h_a/1e3)
    h_p_s.append(h_p/1e3)
plt.figure(figsize=(9, 5))
plt.plot(times[idx_sep+idx_min:]/60, h_a_s, label="Apoapsis")
plt.plot(times[idx_sep+idx_min:]/60, h_p_s, label="Periapsis")
plt.xlabel("Time [min]")
plt.ylabel("Altitude [km]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
print("Inclinaison varies between %.2f and %.2f deg" % (min(inclinaisons), max(inclinaisons)))
print("Apoapsis varies between %.2f and %.2f km" % (min(h_a_s), max(h_a_s)))
print("Periapsis varies between %.2f and %.2f km" % (min(h_p_s), max(h_p_s)))