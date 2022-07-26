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
from scipy import interpolate
from matplotlib import pyplot as plt
import datetime
import glob
import multiprocessing

# Tudatpy imports
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment, environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.astro import element_conversion, time_conversion
from tudatpy.kernel.interface import spice
from tudatpy import util, plotting
spice.load_standard_kernels()

run_sims = False
analyse_res = True

i = 0

def run_orbital_decay(h_p, h_a):
    calendar_date = datetime.datetime(2031, 2, 28, 12, 0, 0)
    julian_date = time_conversion.calendar_date_to_julian_day(calendar_date)
    simulation_start_epoch = time_conversion.julian_day_to_seconds_since_epoch(julian_date)
    simulation_start_epoch = 0

    ## Setup bodies
    current_name = "MAV stage 2"
    bodies_to_create = ["Mars", "Sun"]
    body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Mars", "J2000")
    # Add MCD atmosphere
    from atmosphere import MCD
    mcd = MCD.mcd_interface()
    get_density = mcd.density
    body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_four_dimensional_constant_temperature(
        get_density,
        constant_temperature=215,
        specific_gas_constant=197,
        ratio_of_specific_heats=1.3
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    central_bodies = ["Mars"]
    bodies.create_empty_body(current_name)
    current_body = bodies.get(current_name)
    # Add mass and CD interface
    current_body.set_constant_mass(60)
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
    environment_setup.add_aerodynamic_coefficient_interface(bodies, current_name, coefficient_settings)
    # Add solar radiation pressure interface
    reference_area_radiation = 0.57*2.8*0.25
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Mars"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
        "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
    environment_setup.add_radiation_pressure_interface(bodies, current_name, radiation_pressure_settings)
    bodies_to_propagate = [current_name]

    ## Setup accelerations
    accelerations_on_vehicle = {
        "Mars": [
            propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
            propagation_setup.acceleration.aerodynamic()
        ],
        "Sun": [
            propagation_setup.acceleration.point_mass_gravity(),
            propagation_setup.acceleration.cannonball_radiation_pressure()
        ]
    }
    acceleration_dict = {current_name: accelerations_on_vehicle}
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_dict, bodies_to_propagate, central_bodies
    )

    ## Create initial state
    R_M = spice.get_average_radius("Mars")
    a = (h_p + h_a + 2 * R_M)/2
    r_a, r_p = h_a + R_M, h_p + R_M
    e = 1 - (2 / (r_a/r_p + 1))
    i = np.deg2rad(19.5)
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(a, e, i, 0, 0, 0, spice.get_body_gravitational_parameter("Mars"))

    ## Setup dependent variables
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.apoapsis_altitude(current_name, "Mars"),
        propagation_setup.dependent_variable.periapsis_altitude(current_name, "Mars"),
        propagation_setup.dependent_variable.keplerian_state(current_name, "Mars")
    ]

    ## Setup termination settings
    termination_settings_altitude = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.apoapsis_altitude(current_name, "Mars"),
        limit_value=300e3,
        use_as_lower_limit=True
    )
    def print_t(t):
        global i
        if i % 1000 == 0:
            print("Time %.2f days" % ((t-simulation_start_epoch)/3600/24), end="\r")
        i += 1
        return False
    fake_term_settings = propagation_setup.propagator.custom_termination(print_t)
    termination_max_time_settings = propagation_setup.propagator.time_termination(simulation_start_epoch + 3600*24*365*3)
    termination_settings = propagation_setup.propagator.hybrid_termination([termination_max_time_settings, termination_settings_altitude, fake_term_settings], fulfill_single_condition=True)

    ## Setup propagator settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        propagation_setup.propagator.cowell,
        output_variables=dependent_variables_to_save
    )

    ## Setup integrator settings
    initial_time_step = 15
    minimum_time_step = 3e-5*1000
    maximum_time_step = 600
    tolerance = 1e-7*1000
    coefficients = propagation_setup.integrator.rkf_45
    integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
        simulation_start_epoch,
        initial_time_step,
        coefficients,
        minimum_time_step,
        maximum_time_step,
        relative_error_tolerance=tolerance,
        absolute_error_tolerance=tolerance,
        maximum_factor_increase=1.01,
        minimum_factor_increase=0.01,
        throw_exception_if_minimum_step_exceeded=False,
        assess_termination_on_minor_steps=False
    )

    ## Run simulation
    dynamics_simulator = numerical_simulation.SingleArcSimulator(
        bodies,
        integrator_settings,
        propagator_settings,
        print_dependent_variable_data=False,
        print_state_data=False
    )

    ## Post-processing
    raw_states, raw_dep_vars = dynamics_simulator.state_history, dynamics_simulator.dependent_variable_history
    states = util.result2array(raw_states)
    dep_vars = util.result2array(raw_dep_vars)
    times = (dep_vars[:,0]-simulation_start_epoch)/3600/24
    states = states[:,1:]
    dep_vars = dep_vars[:,1:]
    h_a_s = dep_vars[:,0]/1e3
    h_p_s = dep_vars[:,1]/1e3
    inclinations = np.rad2deg(dep_vars[:,4])
    print()
    print("Re-entered atmosphere after %.2f days"% times[-1])

    ## Save results
    np.savez(sys.path[0]+"/data/optimisation/SA/h_history_%i_%i.npz"%(h_p/1e3, h_a/1e3), times=times, h_a_s=h_a_s, h_p_s=h_p_s, inclinations=inclinations)

if __name__ == "__main__":
    
    if run_sims:
        inputs = []
        altitudes = [305, 315, 325, 335, 345, 355, 365, 375, 385]
        for h_p in altitudes:
            h_p *= 1e3
            for h_a in altitudes:
                h_a *= 1e3
                # Check that apoapsis is higher
                if h_p > h_a:
                    continue
                # Check that results file doesn't exist
                if os.path.isfile(sys.path[0]+"/data/optimisation/SA/h_history_%i_%i.npz"%(h_p/1e3, h_a/1e3)):
                    continue
                print(h_p/1e3, h_a/1e3)
                inputs.append((h_p, h_a))
        input("Press ENTER to run %i inputs..."%len(inputs))
        with multiprocessing.get_context("spawn").Pool(processes=45) as pool:
            outputs = pool.starmap(run_orbital_decay, inputs)

    if analyse_res:
        # Load all files ending in .npz from results folder
        files = glob.glob(sys.path[0]+"/data/optimisation/SA/h_history_*.npz")
        h_a_s, h_p_s, decays = [], [], []
        for file in files:
            results = np.load(file)
            times = results["times"]
            inclinations = results["inclinations"]
            plt.figure(figsize=(9,5))
            plt.plot(times, inclinations)
            idx_0 = np.argwhere(np.fabs((inclinations - inclinations[0])) < 0.1)
            times_repeat = times[idx_0].flatten()
            # Remove values from times_repeat that have less than 0.25 between them
            times_repeat = times_repeat[np.insert(np.diff(times_repeat) > 0.25, 0, False)]

            ylims = plt.ylim()
            plt.vlines(times_repeat, -900, 900, linestyles="dashed", colors="black")
            plt.ylim(ylims)
            plt.xlabel("Time [days]")
            plt.ylabel("Inclination [deg]")
            plt.grid()
            plt.tight_layout()
            plt.savefig(sys.path[0]+"/plots/optimisation/SA/inclination_%s.pdf"%"_".join(file.split("/")[-1].split(".")[0].split("_")[-2:]))
            plt.close()

            plt.figure(figsize=(9,5))
            plt.plot(times, results["h_a_s"], label="Apoapsis")
            plt.plot(times, results["h_p_s"], label="Periapsis")
            plt.xlabel("Time [days]")
            plt.ylabel("Altitude [km]")
            plt.legend(loc="upper right")
            plt.hlines(300, -10, 9999, linestyles="dashed", colors="black")
            plt.xlim(-0.5, times[-1]+0.5)
            plt.grid()
            plt.tight_layout()
            plt.savefig(sys.path[0]+"/plots/optimisation/SA/h_history_%s.pdf"%"_".join(file.split("/")[-1].split(".")[0].split("_")[-2:]))
            plt.close()
            
            _d = file.split(".")[0].split("_")
            h_p, h_a = int(_d[-2]), int(_d[-1])
            print("h_p = %i, h_a = %i, decay after %.2f days"%(h_p, h_a, times[-1]))
            print(" -> , %i repeats at %.2f deg, after %s days"%(len(times_repeat), inclinations[0], str(times_repeat)))
            h_a_s.append(h_a), h_p_s.append(h_p), decays.append(times[-1])

        for x_data, x_label, fname in zip([h_a_s, h_p_s], ["h_p [km]", "h_a [km]"], ["decay_time_vs_h_p", "decay_time_vs_h_a"]):
            plt.figure(figsize=(9,5))
            plt.scatter(x_data, decays)
            plt.xlabel(x_label)
            plt.ylabel("Decay time [days]")
            plt.grid()
            plt.tight_layout()
            plt.savefig(sys.path[0]+"/plots/optimisation/SA/%s.pdf"%fname)
            plt.close()