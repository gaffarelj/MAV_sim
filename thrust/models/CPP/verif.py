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

from thrust.models.tubular import tubular_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM
from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.anchor import anchor_SRM
from thrust.models.spherical import spherical_SRM
from thrust.solid_thrust_multi_stage import SRM_thrust_rk4 as SRM_thrust


SRM_models = []
SRM_models.append(tubular_SRM(R_o=0.24, R_i=0.175, L=1.05))
SRM_models.append(rod_and_tube_SRM(R_o=0.24, R_mid=0.21, R_i=0.175, L=1.05))
SRM_models.append(multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05))
SRM_models.append(anchor_SRM(R_o=1, R_i=0.25, N_a=3, w=0.2, r_f=0.05, delta_s=0.15, L=3))
SRM_models.append(spherical_SRM(R_o=0.175, R_i=0.095))

dt = 1.156191

for SRM_model in SRM_models:
    SRM_thrust_model_1 = SRM_thrust(SRM_model, A_t=0.065, epsilon=45)

    mags = []
    for use_cpp in [False, True]:
        times, magnitudes, b_s, _, M_p_s = SRM_thrust_model_1.simulate_full_burn(dt=dt, use_cpp=use_cpp)
        # Remove negative magnitudes
        magnitudes = [m for m in magnitudes if m > 0]
        # print(magnitudes)
        # print(sum(magnitudes))
        mags.append(sum(magnitudes))
    print("Difference in magnitude: %.3e N (%.2f - %.2f)" % (mags[0] - mags[1], mags[0], mags[1]))