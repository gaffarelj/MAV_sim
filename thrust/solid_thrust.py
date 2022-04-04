import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

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
        T_c=3645,       # [K] chamber temperature (https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/8742205)
        p_a=650*0.4,    # [Pa] ambient pressure (650Pa = "Sea level" on Mars)
        a=5.8e-5,  # [m/s/MPa^n] burning rate coefficient (to use with pressure in MPa) (=0.004202/10**(6*0.31))
        n=0.31,         # [-] burning rate exponent (to use with pressure in MPa)
        rho_p=1854.5,   # [kg/m3] propellant density
        M=0.02414,      # [kg/mol] propellant molar mass FROM CEA # Air Launch versus Ground Launch: a Multidisciplinary Design Optimization Study of Expendable Launch Vehicles on Cost and Performance
        gamma=1.125,    # [-] specific heat ratio of the propellant FROM CEA
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
        # Molar mass of HTPB: 0.0028 kg/mol (https://www.crayvalley.com/docs/default-source/tds/poly-bd-r-45htlo.pdf?sfvrsn=3cdb46f5_6)
        # Molar mass of Al: 0.02698 kg/mol
        # Molar mass of Ammonium perchlorate: 0.11749 kg/mol
        # 70% AP (NH4CLO4) + 16% Al + 14% HTPB
        self.M = M          # = 0.16*0.02698+0.7*0.11749+0.14*0.0028
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
        self.S = self.geometry_model.burning_S(0)          # [m2] burning surface
        self.m_dot = None
        self.last_t = None          # [s] last function call time
        self.saved_burn_times = []
        self.saved_magnitudes = []
        self.saved_m_dot_s = []
        self.magnitude_interpolator = None
        self.m_dot_interpolator = None
        self.p_ratio = None

    def check_At(self, r_max, T_max):
        r_e = np.sqrt(self.A_e/np.pi)
        too_high = r_e > r_max
        T_t = 2 * self.T_c / (self.gamma + 1)
        too_low = T_t > T_max

    def solve_p_c(self, p_c):
        # Solve equation for the chamber pressure
        rho_c = (p_c*self.M/self.R_A/self.T_c)  # [kg/m3] chamber density
        return (self.c_real * (self.rho_p - rho_c) * self.a * self.S / self.A_t)**(1/(1-self.n)) - p_c

    def solve_p_e(self, p_e):
        # Solve equation for the exhaust pressure
        return self.Gamma / np.sqrt(2*self.gamma/(self.gamma-1) * (p_e/self.p_c)**(2/self.gamma) * (1 - (p_e/self.p_c)**((self.gamma-1)/self.gamma))) - self.A_e/self.A_t

    def compute_magnitude(self, time):
        # Compute the current timestep
        if self.last_t is None:
            dt = 1e-12
        else:
            dt = time - self.last_t
        self.last_t = time

        # old_b = self.b
        self.b = self.b + self.r * dt
        # db = self.b - old_b

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
            # self.p_c = fsolve(self.solve_p_c, 1e5)[0]                               # [Pa] combustion chamber pressure
            self.p_c = (self.c_real * self.rho_p * self.a * self.S / self.A_t)**(1/(1-self.n))
            if self.p_ratio is None:
                self.p_e = fsolve(self.solve_p_e, self.p_c/1000)[0]                     # [Pa] exhaust pressure
                self.p_ratio = self.p_e/self.p_c
            else:
                self.p_e = self.p_ratio*self.p_c
            self.V_e = np.sqrt(2*self.gamma/(self.gamma-1) * self.R_A/self.M * self.T_c * (1-(self.p_e/self.p_c)**((self.gamma-1)/self.gamma)))   # [m/s] exhaust velocity
            self.F_T = self.m_dot * self.V_e + (self.p_e - self.p_a) * self.A_e     # [N] thrust
            self.F_T *= self.eta_F_T

            # Update the regression rate
            self.r = self.a * (self.p_c) ** self.n    # [m/s]

            # Update the propellant mass
            self.M_p -= self.m_dot * dt

        # Compute the specific impulse [s]
        self.I_sp = self.V_e / 9.80665
        self.I_sp *= self.eta_Isp

        # Return the thrust
        return self.F_T

    def simulate_full_burn(self, dt=0.01, compute_dep_vars=False, make_interplators=True, filename=None):
        b_s, p_c_s, M_p_s = [], [], []
        if filename is None:
            time = 0
            # Keep computing thrust until the magnitude settles to 0
            while time == 0 or np.sum(self.saved_magnitudes[-2:]) != 0:
                F_T = self.compute_magnitude(time)
                self.saved_magnitudes.append(F_T)
                if compute_dep_vars:
                    b_s.append(self.b)
                    p_c_s.append(self.p_c)
                    #self.saved_Isp_s.append(self.I_sp)
                if compute_dep_vars or make_interplators:
                    self.saved_m_dot_s.append(self.m_dot)
                M_p_s.append(self.M_p)
                self.saved_burn_times.append(time)
                time += dt

            if make_interplators:
                self.magnitude_interpolator = interp1d(self.saved_burn_times, self.saved_magnitudes)
                self.m_dot_interpolator = interp1d(self.saved_burn_times, self.saved_m_dot_s)

            return self.saved_burn_times, self.saved_magnitudes, b_s, p_c_s, M_p_s
        else:
            thrust_results = np.load(filename)
            self.saved_burn_times, self.saved_magnitudes, masses = thrust_results["times"], thrust_results["magnitudes"], thrust_results["masses"]
            mass_diff = np.asarray(np.diff(masses))
            time_diff = np.asarray(np.diff(self.saved_burn_times))
            self.saved_m_dot_s = mass_diff/time_diff
            self.saved_m_dot_s.append(self.saved_m_dot_s[-1])
            
            self.magnitude_interpolator = interp1d(self.saved_burn_times, self.saved_magnitudes)
            self.m_dot_interpolator = interp1d(self.saved_burn_times, self.saved_m_dot_s)

            return self.saved_burn_times, self.saved_magnitudes, b_s, p_c_s, masses