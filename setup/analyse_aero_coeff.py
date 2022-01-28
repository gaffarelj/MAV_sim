import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize, interpolate

# Atmospheric conditions
hs =   [100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]
rhos = [1.82031e-07, 4.29462e-09, 1.94534e-10, 2.15307e-11, 3.60822e-12, 8.64996e-13, 2.74321e-13, 1.01443e-13, 4.26962e-14, 1.10285e-14, 4.70352e-15, 2.97560e-15, 2.29431e-15]
ps =   [4.58579e-03, 1.07168e-04, 7.27184e-06, 9.83263e-07, 2.05882e-07, 6.32226e-08, 2.59061e-08, 1.29177e-08, 7.73174e-09, 4.31813e-09, 3.26268e-09, 2.72788e-09, 2.35800e-09]
Ts =   [128.643,     133.166,     172.404,     178.405,     179.196,     179.330,     179.369,     179.380,     179.382,     179.383,     179.383,     179.382,     179.383]
Vs =   [3503.34,     3490.86,     3478.51,     3466.29,     3454.20,     3442.23,     3430.39,     3418.67,     3407.07,     3384.21,     3361.81,     3339.85,     3318.31]
fracs = [
    np.array([93.270, 2.445, 2.513, 0.965, 0.607, 0.200])/100,
    np.array([81.674, 6.398, 4.972, 3.851, 2.509, 0.596])/100,
    np.array([65.826, 11.923, 4.910, 7.796, 8.533, 1.012])/100,
    np.array([44.654, 16.782, 3.754, 10.801, 22.822, 1.186])/100,
    np.array([24.723, 17.401, 2.297, 11.148, 43.375, 1.055])/100,
    np.array([11.796, 14.477, 1.210, 9.397, 62.345, 0.775])/100,
    np.array([5.102, 10.644, 0.582, 7.033, 76.132, 0.507])/100,
    np.array([1.900, 7.190, 0.249, 4.862, 85.498, 0.301])/100,
    np.array([0.648, 4.581, 0.098, 3.174, 91.333, 0.166])/100,
    np.array([0.074, 1.758, 0.015, 1.271, 96.835, 0.047])/100,
    np.array([0.012, 0.679, 0.003, 0.501, 98.792, 0.014])/100,
    np.array([0.004, 0.280, 0.001, 0.203, 99.507, 0.005])/100,
    np.array([0.003, 0.132, 0.001, 0.088, 99.775, 0.002])/100
]
Ms = [19.5229, 18.7224, 15.7385, 14.3184, 12.9100, 11.7434, 10.9463, 10.4315, 10.1257, 9.8292, 9.6898, 9.6002, 9.5286]

# Drag coefficients
CDs = {
    100: 1.59860,
    125: 1.77670,
    150: 1.75737,
    175: 1.74900,
    200: 1.74924,
    225: 1.70594,
    250: 1.69850,
    275: 1.71255,
    300: 1.68698,
    350: 1.67177,
    400: 1.67533,
    450: 1.68058,
    500: 1.67640
}

def func(x, a, b, c):
    return a * np.exp(-x / b) + c

import time
t0 = time.time()

popt, pcov = optimize.curve_fit(func, list(CDs.keys()), list(CDs.values()), p0=[1.2, 70, 1.65])

print(time.time() - t0)
print(popt)
print(pcov)
altitudes = np.arange(75, 600, 0.1)
CDs_fit = [func(h, *popt) for h in altitudes]

cs = interpolate.interp1d(list(CDs.keys()), list(CDs.values()), kind="quadratic", bounds_error=False, fill_value="extrapolate")
CDs_interp = [cs(h) for h in altitudes]

plt.scatter(CDs.keys(), CDs.values(), label="DSMC coefficients")
plt.plot(altitudes, CDs_interp, color="orange", label="Spline interpolation")
plt.plot(altitudes, CDs_fit, color="green", label="Exponential fit: %.2f * np.exp(-x / %.2f) + %.2f" % tuple(popt))
plt.xlabel("Altitude [km]"), plt.ylabel("Drag coefficient [-]")
plt.grid()
plt.legend()
plt.show()


# Also plot atmospheric properties as function of altitude + coefficients as function of Mach number

# plt.scatter(hs, rhos)
# plt.xlabel("Altitude [km]"), plt.ylabel("Atmospheric density [kg/m$^3$]")
# plt.yscale("log")
# plt.grid()
# plt.show()

# plt.scatter(hs, ps)
# plt.xlabel("Altitude [km]"), plt.ylabel("Atmospheric pressure [Pa]")
# plt.yscale("log")
# plt.grid()
# plt.show()

# plt.scatter(hs, Ts)
# plt.xlabel("Altitude [km]"), plt.ylabel("Temperature [K]")
# plt.grid()
# plt.show()

# plt.scatter(hs, Vs)
# plt.xlabel("Altitude [km]"), plt.ylabel("Velocity [m/s]")
# plt.grid()
# plt.show()

# plt.scatter(hs, Ms)
# plt.xlabel("Altitude [km]"), plt.ylabel("Mach number [-]")
# plt.grid()
# plt.show()

plt.scatter(CDs.values(), Ms)
plt.xlabel("Drag coefficient [-]"), plt.ylabel("Mach number [-]")
plt.grid()
plt.show()
