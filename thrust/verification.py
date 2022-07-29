import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
    
from thrust.solid_thrust import SRM_thrust
from thrust.models.tubular import tubular_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.multi_fin import multi_fin_SRM

# https://github.com/reilleya/openMotor/blob/staging/test/data/regression/simple/motor.ric

def load_data(fname):
    path = sys.path[0] + "/data/openMotor_results/%s.csv" % fname
    df = pd.read_csv(path)
    required_columns = ['Time(s)', 'Chamber Pressure(MPa)', 'Thrust(N)', 'Propellant Mass(G1;kg)', 'Mass Flow(G1;kg/s)', 'Regression Depth(G1;m)', 'Nozzle Exit Pressure(MPa)']
    missing_columns = set(required_columns).difference(set(df.columns))
    if len(missing_columns) != 0:
        raise ValueError("The openMotor result file is missing some columns:", missing_columns)
    return [df[col].values for col in required_columns]

tubular_GEO = tubular_SRM(
    R_o=0.04153,
    R_i=0.01588,
    L=0.1397
)
SRM_thrust_model_tubular = SRM_thrust(
    tubular_GEO,
    A_t=0.00015327962,
    epsilon=6.25,
    a=3.5952e-05,
    n=0.3273,
    T_c=2800,
    p_a=101325,
    rho_p=1670,
    M=0.02367,
    gamma=1.21,
    eta_F_T=0.9
)

rod_and_tube_GEO = rod_and_tube_SRM(
    R_o=0.1841503683007366/2,
    R_mid=0.10731521463042927/2,
    R_i=0.07683515367030735/2,
    L=0.9144018288036577
)
SRM_thrust_model_rod_and_tube = SRM_thrust(
    rod_and_tube_GEO,
    A_t=0.00137951662,
    epsilon=8.28742,
    a=3.517e-05,
    n=0.3273,
    T_c=3500,
    p_a=101325,
    rho_p=1680,
    M=0.02367,
    gamma=1.21,
    eta_F_T=0.9
)

# multi_fin_GEO = multi_fin_SRM(
#     R_o=1.0/2,
#     R_i=0.65/2,
#     N_f=8,
#     w_f=0.15/2,
#     L_f=0.3/2,
# L=3.0
# )
# SRM_thrust_model_multi_fin = SRM_thrust(
#     multi_fin_GEO,
#     A_t=0.065,
#     epsilon=12.083,
#     a=0.004177,
#     n=0.059,
#     T_c=1520,
#     p_a=101325,
#     rho_p=1750,
#     M=0.0399,
#     gamma=1.1361,
#     eta_F_T=0.9
# )

SRM_models = [SRM_thrust_model_tubular, SRM_thrust_model_rod_and_tube]
model_names = ["tubular", "rod_and_tube"]

for model_name, SRM_model in zip(model_names, SRM_models):
    # Plot the geometry of the SRM model
    SRM_model.geometry_model.plot_geometry()

    # Simulate the full SRM burn
    times, magnitudes, b_s, p_c_s, M_p_s = SRM_model.simulate_full_burn(dt=1e-3, compute_dep_vars=True)
    p_c_s = np.asarray(p_c_s)/1e6
    times_OM, p_c_OM, magnitudes_OM, mass_OM, mdot_OM, b_OM, p_e_OM = load_data(model_name)

    # Compile validation and simulation data
    names = ["Thrust [N]", "Propellant mass [kg]"]#, "Regression depth [m]", "Chamber pressure [MPa]"]
    magnitudes = [0 if magnitude < 0 else magnitude for magnitude in magnitudes]
    p_c_s[-1] = 0
    values = [magnitudes, M_p_s]#, b_s, p_c_s]
    values_OM = [magnitudes_OM, mass_OM]#, b_OM, p_c_OM]

    # Plot
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,8))
    fig, ((ax1), (ax2)) = plt.subplots(1, 2, figsize=(9,3.5))
    # axes = [ax1, ax2, ax3, ax4]
    axes = [ax1, ax2]
    # fig.suptitle("%s geometry" % model_name.replace("_", " ").capitalize())
    for name, value, value_OM, ax in zip(names, values, values_OM, axes):
        ax.plot(times, value, label="Simulation", linestyle="dotted")
        ax.plot(times_OM, value_OM, label="openMotor", linestyle="dotted")
        ax.set_xlabel("Time [s]"), ax.set_ylabel(name)
        ax.grid(), ax.legend()
    plt.tight_layout()
    plt.savefig(sys.path[0] + "/plots/%s_geometry_verif_cropped.pdf" % model_name)