import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

use_all_propellants = True

# Data taken from the TRP reader, p.494 (p.500 in pdf)
# CTBP14/Al16/AP70
name_s =    ["STAR 13B", "STAR 17", "STAR 24", "STAR 17A", "STAR 24C"]
thrust_s =  [7600, 10900, 18500, 16000, 20300]
I_sp_s =    [285, 296.2, 282.9, 286.7, 282.3]
p_c_s =     [5.67e6, 5.54e6, 3.35e6, 4.62e6, 3.75e6]
M_p_s =     [41, 70, 200, 112, 220]
M_inrt_s =  [5.81, 9.43, 18.3, 13.4, 19.7]
L_s =       [0.638, 0.687, 1.03, 0.981, 0.64]
D_s =       [0.3454, 0.442, 0.622, 0.442, 0.622]
t_b_s =     [14.8, 17.6, 29.6, 19.4, 29.6]
epsilon_s = [49.8, 60.7, 37.8, 53.2, 37.1]

cols = ["Name", "Thrust [N]", "I_sp [s]", "p_c [Pa]", "m_prop [kg]", "m_innert [kg]", "L [m]", "D [m]", "t_b [s]", "epsilon [-]"]

df_CTBP = pd.DataFrame(np.array([name_s, thrust_s, I_sp_s, p_c_s, M_p_s, M_inrt_s, L_s, D_s, t_b_s, epsilon_s]).T, columns=cols)

# HTBP11/Al18/AP71
name_s =    ["STAR 31", "STAR 75", "STAR 27", "STAR 30E", "STAR 37FM", "STAR 27H", "STAR 30Eb", "STAR 26", "STAR 30C", "STAR 30BP", "STAR 27b", "STAR 26B", "STAR 30Cb", "STAR 48A short", "STAR 48A long", "STAR 48B", "STAR 48V"]
thrust_s =  [82.3e3, 200.17e3, 27e3, 35.2e3, 47.85e3, 20.68e3, 35.1e3, 33.4e3, 32.5e3, 26.6e3, 25.4e3, 34.62e3, 32.47e3, 77.18e3, 78.96e3, 67.2e3, 68.6e3]
I_sp_s =    [293.5, 290, 272.4, 292.8, 291.8, 291.4, 290.4, 271, 291.8, 293, 287.9, 272.4, 288.8, 285.3, 291.9, 286, 292.1]
p_c_s =     [49.1e5, 42.47e5, 38.82e5, 36.9e5, 37.23e5, 41.09e5, 37e5, 39.6e5, 38.1e5, 35.4e5, 38.8e5, 42.95e5, 38.06e5, 37.44e5, 37.44e5, 39.9e5, 39.9e5]
M_p_s =     [1286, 7503, 334, 631, 1057, 337.8, 631, 231, 591, 515, 334, 238, 590.8, 2430, 2430, 2010, 2010]
M_inrt_s =  [108, 562.9, 27.49, 42.5, 91.62, 29.98, 42.5, 38.8, 41.1, 37.7, 27.5, 23.4, 38.46, 134.6, 151.5, 124, 154]
L_s =       [2.87, 2.59, 1.303, 1.68, 1.68, 1.219, 1.69, 0.838, 1.633, 1.829, 1.24, 0.841, 1.49, 2.032, 2.235, 2.032, 2.075]
D_s =       [0.764, 1.905, 0.6934, 0.762, 0.933, 0.6934, 0.762, 0.66, 0.7612, 0.762, 0.6934, 0.663, 0.762, 1.2426, 1.2426, 1.24, 1.24]
t_b_s =     [45, 105, 34.5, 49.3, 63.7, 46.3, 51.1, 17.8, 51, 54, 37.3, 17.8, 51, 87.2, 87.2, 84.1, 84.1]
epsilon_s = [58, 17.7, 48.8, 8.6, 48.2, 81.7, 36.9, 16.7, 63.2, 73.7, 48.8, 17.8, 46.4, 32.2, 43.1, 39.6, 54.8]

df_HTBP = pd.DataFrame(np.array([name_s, thrust_s, I_sp_s, p_c_s, M_p_s, M_inrt_s, L_s, D_s, t_b_s, epsilon_s]).T, columns=cols)

if use_all_propellants:
    df = pd.concat([df_CTBP, df_HTBP], ignore_index=True)
else:
    df = df_CTBP

for col in df.columns:
    t = str if col == "Name" else float
    df[col] = df[col].astype(t)


df["Area [m2]"] = np.pi*df["D [m]"]*df["L [m]"]
print(df)

# plt.scatter(df["Area [m2]"], df["p_c [Pa]"])
# plt.xlabel("Area [m2]"), plt.ylabel("p_c [Pa]")
# plt.grid()
# plt.tight_layout()
# plt.show()

sns.pairplot(df, vars=df.columns[1:], corner=True)#, kind="kde")
plt.show()