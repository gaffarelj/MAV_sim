import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

import numpy as np
from matplotlib import pyplot as plt

#####################################
### Test the tubular SRM geometry ###
#####################################
from models.tubular import tubular_SRM

# Define the geometry
R_o = 1
R_i = 0.25
L = 3
tubular_test = tubular_SRM(R_o, R_i, L)

# Compute the burning area
b_s = np.arange(0, R_o-R_i, 0.0001)
P_s = [tubular_test.burning_S(b) for b in b_s]

# Plot the initial geometry
tubular_test.plot_geometry()

# Save the geometry burn over time
for b in np.linspace(0, R_o-R_i-0.0001, 50):
    tubular_test.plot_geometry(b, save=sys.path[0]+"/thrust/burn_visu/tubular_%.4f.png" % b)

# Plot the burning area over burnt thickness
plt.plot(b_s, P_s)
plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning area [m$^2$]")
plt.grid(), plt.tight_layout()
plt.show()

##########################################
### Test the rod and tube SRM geometry ###
##########################################
from models.rod_and_tube import rod_and_tube_SRM

# Define the geometry
R_o = 1
R_mid = 0.75
R_i = 0.25
L = 3
rod_and_tube_test = rod_and_tube_SRM(R_o, R_mid, R_i, L)

# Compute burning perimeter
b_s = np.arange(0, max(R_o - R_mid, R_i), 0.0001)
P_s = [rod_and_tube_test.burning_S(b) for b in b_s]

# Plot the initial geometry
rod_and_tube_test.plot_geometry()

# Save the geometry burn over time
for b in np.linspace(0, max(R_o - R_mid, R_i)-0.0001, 50):
    rod_and_tube_test.plot_geometry(b, save=sys.path[0]+"/thrust/burn_visu/rod_and_tube_%.4f.png" % b)

# # Plot the burning area over burnt thickness
# plt.plot(b_s, P_s)
# plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning area [m$^2$]")
# plt.grid(), plt.tight_layout()
# plt.show()

##########################################
### Test the multi-fin SRM geometry ###
##########################################
from models.multi_fin import multi_fin_SRM

# Define the geometry
R_o = 1
R_i = 0.75
N_f = 10
w_f = 0.15
L_f = 0.35
L = 3
multi_fin_test = multi_fin_SRM(R_o, R_i, N_f, w_f, L_f, L)

# Compute burning perimeter
b_s = np.arange(0, R_o-R_i, 0.0001)
P_s = [multi_fin_test.burning_S(b) for b in b_s]

# Plot the initial geometry
multi_fin_test.plot_geometry()

# Save the geometry burn over time
for b in np.linspace(0, R_o-R_i-0.0001, 50):
    multi_fin_test.plot_geometry(b, save=sys.path[0]+"/thrust/burn_visu/multi_fin_%.4f.png" % b)

# Plot the burning area over burnt thickness
plt.plot(b_s, P_s)
plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning area [m$^2$]")
plt.grid(), plt.tight_layout()
plt.show()