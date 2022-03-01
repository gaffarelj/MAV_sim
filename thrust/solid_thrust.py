import numpy as np
from scipy.optimize import fsolve

import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

########################################
### PROPELLANT TO BE USED: TP-H-3544 ###
# https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/8742205
########################################

# TODO: find data of remaining elements

class SRM_thrust:

    def __init__(self,
        geometry_model,
        A_t=0.0019635,  # [m2] throat area
        epsilon=50.2,   # [-] area ratio (exhaust/throat)
        T_c=3645,       # [K] combustion temperature (https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/8742205)
        p_a=650,        # [Pa] ambient pressure (650Pa = "Sea level" on Mars)
        a=0.004202,     # [...] burning rate coefficient (to use with pressure in MPa)
        n=0.31,         # [-] burning rate exponent (to use with pressure in MPa)
        rho_p=1854.5,   # [kg/m3] propellant density
        M=0.4785598,    # [kg/mol] propellant molar mass
        gamma=1.175,    # [-] specific heat ratio of the propellant
        eta_Isp=0.95,   # [-] Specific impulse efficiency
        eta_c=0.93,     # [-] Combustion efficiency
        eta_F_T=0.95):  # [-] Thrust efficiency

        self.geometry_model = geometry_model
        self.A_t = A_t
        self.A_e = A_t * epsilon
        self.T_c = T_c              
        self.p_a = p_a

        # Propellant data for CTPB TP-H-3062 (CTPB14/Al16/AP70), similar to TP-H-3544 (TP-H-3544 is HTPB [does not self-degrade, not even in vacuum of space])
        # REF: p.12 of http://31.186.81.235:8080/api/files/view/38619.pdf
        self.a = a              
        self.n = n              
        # Molar mass of HTPB: 2.8 kg/mol (https://www.crayvalley.com/docs/default-source/tds/poly-bd-r-45htlo.pdf?sfvrsn=3cdb46f5_6)
        # Molar mass of Al: 0.02698 kg/mol
        # Molar mass of Ammonium perchlorate: 0.11749 kg/mol
        self.M = M          # = 0.16*0.02698+0.7*0.11749+0.14*2.8
        self.R_A = 8.314    # [J/mol/K] gas constant
        # 1749 kg/m3 for TP-H-3062 in DAMNOGLYOXIME AND DAMNOFURAZAN IN PROPELLANTS BASED ON AMMONUMPERCHLORATE
        # Mars Ascent Vehicle—Propellant Aging: TP-H-3544 has higher density than TP-H-3062
        # Rocket propulsion elements (p.479): 1854.5 kg/m3
        self.rho_p = rho_p
        self.gamma = gamma

        #########
        #### Rocket Propulsion Elements: Figure 15-6: adding 0-20% of ultra-fine aluminum increases regression rate by 0-75%
        #########

        # Efficiencies
        self.eta_Isp = eta_Isp      # TRP Table 2 p.271
        self.eta_F_T = eta_F_T
        self.eta_c = eta_c

        self.M_p = self.geometry_model.get_V_p() * self.rho_p   # [kg] total propellant mass
        self.M_innert = 0.2833*self.M_p**0.789+10               # [kg] innert mass with titanium casing (with R^2=0.985) # REF: p.497 TRP reader # added 10 kg for TVC, based on results from Mars_Ascent_Vehicle_MAV_Propulsion_Subsystems_Design

        # Compute invariables
        self.Gamma = np.sqrt(self.gamma) * (2/(self.gamma+1))**((self.gamma+1)/(2*(self.gamma-1)))  # [-] Vandenkerckhove function
        self.c_ideal = np.sqrt(self.R_A/self.M*self.T_c)/self.Gamma                                 # [m/s] ideal characteristic velocity
        self.c_real = self.eta_c * self.c_ideal                                                     # [m/s] real characteristic velocity

        # Initial assumptions
        self.r = 1e-9       # [m/s] initial regression rate, close to 0
        self.b = 0          # [m] burnt thickness
        self.S = 0          # [m2] burning surface

        self.m_dot = None
        self.last_t = None          # [s] last function call time

    def solve_p_c(self, p_c):
        # Solve equation for the chamber pressure
        rho_c = (p_c*self.M/self.R_A/self.T_c)  # [kg/m3] chamber density
        return (self.c_real * (self.rho_p - rho_c) * self.a * self.S / self.A_t / 1e6)**(1/(1-self.n)) * 1e6 - p_c

    def solve_p_e(self, p_e):
        # Solve equation for the exhaust pressure
        return self.Gamma / np.sqrt(2*self.gamma/(self.gamma-1) * (p_e/self.p_c)**(2/self.gamma) * (1 - (p_e/self.p_c)**((self.gamma-1)/self.gamma))) - self.A_e/self.A_t

    def compute_magnitude(self, time):
        # Compute the current timestep
        if self.last_t is None:
            dt = 1e-6
        else:
            dt = time - self.last_t
        self.last_t = time
        self.b = self.b + self.r * dt

        # Get the current burning surface
        try:
            self.S = self.geometry_model.burning_S(self.b)
        # Return 0 if the burn thickness is too high
        except ValueError:
            return 0

        # Stop the thrust if there is no surface burning anymore
        if self.S == 0:
            self.m_dot, self.p_c, self.p_e, self.V_e, self.F_T, self.r = 0, 0, 0, 0, 0, 0
        else:
            # Compute propulsion characteristics at this time step
            self.m_dot = self.S * self.r * self.rho_p                               # [kg/s] mass flow
            self.p_c = fsolve(self.solve_p_c, 1e5)[0]                               # [Pa] combustion chamber pressure
            self.p_e = fsolve(self.solve_p_e, self.p_c/1000)[0]                     # [Pa] exhaust pressure
            self.V_e = np.sqrt(2*self.gamma/(self.gamma-1) * self.R_A/self.M * self.T_c * (1-(self.p_e/self.p_c)**((self.gamma-1)/self.gamma)))   # [m/s] exhaust velocity
            self.F_T = self.m_dot * self.V_e + (self.p_e - self.p_a) * self.A_e     # [N] thrust
            self.F_T *= self.eta_F_T

            # Update the regression rate
            self.r = self.a * (self.p_c/1e6) ** self.n    # [m/s]

            # Update the propellant mass
            self.M_p -= self.m_dot * dt

        # Compute the specific impulse [s]
        self.I_sp = self.V_e / 9.80665
        self.I_sp *= self.eta_Isp

        # Return the thrust
        return self.F_T

    def get_Isp(self, time=None):
        return self.I_sp

    def get_m_dot(self, time):
        if self.m_dot is None or self.F_T == 0:
            return 0
        return -np.fabs(self.m_dot)