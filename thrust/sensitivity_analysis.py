import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Classic imports
import numpy as np
from matplotlib import pyplot as plt

# Custom imports
from thrust.solid_thrust import SRM_thrust
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.spherical import spherical_SRM

# Define SRM input names
model_inputs_names = [
    "$A_t$ [m$^2$]",
    "$\epsilon$ [-]",
    "$T_c$ [K]",
    "$p_a$ [Pa]",
    "$a$ [m/s/MPa]",
    "$n$ [-]",
    "$\\rho_p$ [kg/m$^3$]",
    "$M$ [kg/mol]",
    "$\gamma$ [-]",
    "$eta_{I_{sp}}$ [-]",
    "$eta_c$ [-]",
    "$eta_{F_T} [-]$"
]

# Define SRM mean inputs for on-at-a-time analysis
mean_model_inputs = [0.0019635, 50.2, 3645, 650, 0.004202, 0.31, 1854.5, 0.4785598, 1.175, 0.95, 0.93, 0.95]

# Define multiples to apply to each input individually
model_inputs_multiples = [
    2**np.linspace(-4.5, 3.5, 17),      # throat area
    2**np.linspace(-1, 1, 13),          # area ratio
    2**np.linspace(-0.25, 0.25, 9),     # chamber temperature
    [1],                                # ambient pressure
    2**np.linspace(-0.15, 0.15, 9),     # burning rate coefficient
    2**np.linspace(-0.15, 0.15, 9),     # burning rate exponent
    2**np.linspace(-0.25, 0.25, 13),    # propellant density
    2**np.linspace(-0.15, 0.15, 9),     # propellant molar mass
    2**np.linspace(-0.2, 0.2, 13),      # specific heat ratio of the propellant
    [1],                                # specific impulse efficiency
    [1],                                # combustion efficiency
    [1]                                 # thrust efficiency
]

# Compute all model inputs
model_inputs = [mean_model_inputs[i]*np.array(model_inputs_multiples[i]) for i in range(len(mean_model_inputs))]

# Define thrust model
# SRM_geometry = multi_fin_SRM(R_o=0.275, R_i=0.15, N_f=12, w_f=0.0225, L_f=0.075, L=1.125)
# SRM_geometry = spherical_SRM(R_o=0.165, R_i=0.05)
SRM_geometry = rod_and_tube_SRM(R_o=0.28, R_mid=0.14, R_i=0.135, L=1.125)

# Define plotting colors and linestyles
colors = ["C%i"%i for i in range(10)]
linestyles = ["solid", "dotted", "dashed"]
i_c, i_ls = 0, 0

# Loop trough all inputs and plot the resulting thrust profile
for i, model_input_multiples in enumerate(model_inputs_multiples):
    if len(model_input_multiples) > 1:
        fig = plt.figure(figsize=(14,7))
        ax = fig.add_subplot()

        for j, specific_input in enumerate(model_inputs[i]):
            print("Running input %i/%i variation %i/%i..." % (i+1, len(model_inputs_multiples), j+1, len(model_inputs[i])), end="\r")
            current_inputs = mean_model_inputs.copy()
            current_inputs[i] = specific_input
            SRM_thrust_model = SRM_thrust(SRM_geometry, *current_inputs)

            burn_times, magnitudes, b_s, p_c_s, M_p_s = SRM_thrust_model.simulate_full_burn(dt=0.01)

            ax.plot(burn_times, np.asarray(magnitudes)/1e3, label="%s = %.3e" % (model_inputs_names[i], specific_input), color=colors[i_c], linestyle=linestyles[i_ls])

            i_c += 1
            if i_c == len(colors):
                i_c = 0
                i_ls += 1

        print()
        ax.set_xlabel("Burn time [s]"), ax.set_ylabel("Thrust magnitude [kN]")
        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.show()