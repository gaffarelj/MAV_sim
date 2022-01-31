import numpy as np
from matplotlib import pyplot as plt

# Define geometry
R_o = 1
R_i = 0.25

# Check conditions
c_1 = R_i < R_o
print(c_1)

# Compute burning perimeter
b_s = np.arange(0, R_o-R_i, 0.0001)
P_s = 2*np.pi * (R_i + b_s)

# Plot
plt.plot(b_s, P_s, label="Total")

plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning perimeter [m]")
# plt.legend()
plt.grid(), plt.tight_layout()
plt.show()