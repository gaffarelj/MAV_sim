import numpy as np
from scipy.optimize import fsolve


def solve_p_c(p_c):
    # Solve equation for the chamber pressure
    return (c_real * (rho_p - (p_c*M/R_A/T_c) ) * a * S / A_t)**(1/(1-n)) - p_c

def solve_p_e(p_e):
    # Solve equation for the exhaust pressure
    return Gamma / np.sqrt(2*gamma/(gamma-1) * (p_e/p_c)**(2/gamma) * (1 - (p_e/p_c)**((gamma-1)/gamma))) - A_e/A_t


# Given values
gamma = 1.19    # [-] specific heat ratio
a = 0.004       # [-] burning rate coefficient
n = 0.3         # [-] burning rate exponent
rho_p = 1800    # [kg/m3] propellant density
A_t = 2e-4      # [m2] throat area
A_e = 4e-3      # [m2] exhaust area
M = 0.031       # [kg/mol] molar mass
R_A = 8.314     # [J/mol/K] gas constant
eta_c = 0.93    # [-] combustion efficiency
T_c = 3400      # [K] combustion temperature
p_a = 1e5       # [Pa] ambient pressure

dt = 0.001  # [s] time step
S = 2       # [m2] burning surface

# Initial assumptions
dS = 0.5    # [m2] difference of surface burned between t_i and t_i-1
r = 0.001   # [m/s] initial regression rate, close to 0

# Compute invariables
Gamma = np.sqrt(gamma) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))   # [-] Vandenkerckhove function
c_ideal = np.sqrt(R_A/M*T_c)/Gamma                                  # [m/s] ideal characteristic velocity
c_real = eta_c * c_ideal                                            # [m/s] real characteristic velocity

# Compute values at this time step
m_dot = dS * r * rho_p                              # [kg/s] mass flow
p_c = fsolve(solve_p_c, 1e5)[0]                     # [Pa] combustion chamber pressure
p_e = fsolve(solve_p_e, p_c/100)[0]                 # [Pa] exhaust pressure
V_e = np.sqrt(2*gamma/(gamma-1) * R_A/M * T_c * (1-(p_e/p_c)**((gamma-1)/gamma)))   # [m/s] exhaust velocity
F_T = m_dot * V_e + (p_e - p_a) * A_e               # [N] thrust

print(F_T)
