import numpy as np
from matplotlib import pyplot as plt

# Define geometry
R_o = 1
R_i = 0.35
N_f = 5
w_f = 0.1
L_f = 0.25

# Check conditions
c_1 = w_f/2 < R_o - R_i
c_2 = L_f < R_i
c_3 = np.sin(np.pi/N_f) * (w_f + 2 * (R_i - L_f)) >= 0
print(c_1, c_2, c_3)

# Compute burning perimeter
b_s = np.arange(0, R_o-R_i, 0.0001)
P_s = []
for b in b_s:
    P_tube = 2*np.pi * (R_i + b)
    P_fin = 0
    if b < w_f/2:
        P_fin = 2 * N_f * (L_f - b)
    P_s.append(P_tube + P_fin)

# Plot
plt.plot(b_s, P_s, label="Total")

plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning perimeter [m]")
# plt.legend()
plt.grid(), plt.tight_layout()
plt.show()