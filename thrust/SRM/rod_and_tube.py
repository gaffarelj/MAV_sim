import numpy as np
from matplotlib import pyplot as plt

# Define geometry
R_o = 1
R_mid = 0.75
R_i = 0.25

# Check conditions
c_1 = R_i < R_mid
c_2 = R_mid < R_o
print(c_1, c_2)

# Compute burning perimeter
b_s = np.arange(0, max(R_o - R_mid, R_i), 0.0001)
P_s = []
for b in b_s:
    P_inner, P_outer = 0, 0
    if R_mid + b < R_o:
        P_outer = 2 * np.pi * (R_mid + b)
    if R_i - b > 0:
        P_inner = 2 * np.pi * (R_i - b)
    P_s.append(P_inner + P_outer)

# Plot
plt.plot(b_s, P_s, label="Total")

plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning perimeter [m]")
# plt.legend()
plt.grid(), plt.tight_layout()
plt.show()