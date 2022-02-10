import numpy as np
from scipy.optimize import fsolve

import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))


class thrust:

    def __init__(self, geometry_model, rho_p=1800, A_t=0.0019635, A_e=0.1452):
        self.geometry_model = geometry_model
        self.gamma = 1.19   # [-] specific heat ratio
        self.a = 0.004      # [-] burning rate coefficient
        self.n = 0.3        # [-] burning rate exponent
        self.rho_p = rho_p  # [kg/m3] propellant density
        self.A_t = A_t      # [m2] throat area
        self.A_e = A_e      # [m2] exhaust area
        self.M = 0.031      # [kg/mol] molar mass
        self.R_A = 8.314    # [J/mol/K] gas constant
        self.eta_c = 0.93   # [-] combustion efficiency
        self.T_c = 3400     # [K] combustion temperature
        self.p_a = 650      # [Pa] ambient pressure
        self.last_t = None  # [s] last function call time

        # Compute invariables
        self.Gamma = np.sqrt(self.gamma) * (2/(self.gamma+1))**((self.gamma+1)/(2*(self.gamma-1)))  # [-] Vandenkerckhove function
        self.c_ideal = np.sqrt(self.R_A/self.M*self.T_c)/self.Gamma                                 # [m/s] ideal characteristic velocity
        self.c_real = self.eta_c * self.c_ideal                                                     # [m/s] real characteristic velocity

        # Initial assumptions
        self.r = 1e-9       # [m/s] initial regression rate, close to 0
        self.b = 0          # [m] burnt thickness
        self.S = 0          # [m2] burning surface

    def solve_p_c(self, p_c):
        # Solve equation for the chamber pressure
        rho_c = (p_c*self.M/self.R_A/self.T_c)  # [kg/m3] chamber density
        return (self.c_real * (self.rho_p - rho_c) * self.a * self.S / self.A_t / 1e6)**(1/(1-self.n)) * 1e6 - p_c

    def solve_p_e(self, p_e):
        # Solve equation for the exhaust pressure
        return self.Gamma / np.sqrt(2*self.gamma/(self.gamma-1) * (p_e/self.p_c)**(2/self.gamma) * (1 - (p_e/self.p_c)**((self.gamma-1)/self.gamma))) - self.A_e/self.A_t

    def magnitude(self, time):
        # Compute the current timestep
        if self.last_t is None:
            dt = 1e-3
        else:
            dt = time - self.last_t
        self.last_t = time
        self.b = self.b + self.r * dt

        # Get the current burning surface
        try:
            S = self.geometry_model.burning_S(self.b)
        # Return 0 if the burn thickness is too high
        except ValueError:
            return 0

        # Compute the change in burning surface
        dS = max(S - self.S, 0)
        self.S = S

        # Stop the thrust if there is no surface burning anymore
        if self.S == 0:
            self.m_dot, self.p_c, self.p_e, self.V_e, self.F_T, self.r = 0, 0, 0, 0, 0, 0
        else:
            # Compute propulsion characteristics at this time step
            self.m_dot = dS * self.r * self.rho_p                                   # [kg/s] mass flow
            self.p_c = fsolve(self.solve_p_c, 1e5)[0]                               # [Pa] combustion chamber pressure
            self.p_e = fsolve(self.solve_p_e, self.p_c/1000)[0]                     # [Pa] exhaust pressure
            self.V_e = np.sqrt(2*self.gamma/(self.gamma-1) * self.R_A/self.M * self.T_c * (1-(self.p_e/self.p_c)**((self.gamma-1)/self.gamma)))   # [m/s] exhaust velocity
            self.F_T = self.m_dot * self.V_e + (self.p_e - self.p_a) * self.A_e     # [N] thrust

            # Update the regression rate
            self.r = self.a * (self.p_c/1e6) ** self.n    # [m/s]

        # Return the thrust
        return self.F_T