import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

class SRM_thrust_rk4:

    def __init__(self,
        geometry_model,
        A_t=0.0019635,  # [m2] throat area
        epsilon=50.2,   # [-] area ratio (exhaust/throat)
        T_c=3645,       # [K] chamber temperature (https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/8742205)
        p_a=650*0.4,    # [Pa] ambient pressure (650Pa = "Sea level" on Mars)
        a=0.004202/10**(6*0.31), # [m/s/MPa^n] burning rate coefficient (to use with pressure in MPa) (=0.004202/10**(6*0.31))
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

        self.a = a              
        self.n = n              
        self.M = M          # = 0.16*0.02698+0.7*0.11749+0.14*0.0028
        self.R_A = 8.314    # [J/mol/K] gas constant
        self.rho_p = rho_p
        self.gamma = gamma

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

    def solve_p_e(self, p_e):
        # Solve equation for the exhaust pressure
        return self.Gamma / np.sqrt(2*self.gamma/(self.gamma-1) * (p_e/self.p_c)**(2/self.gamma) * (1 - (p_e/self.p_c)**((self.gamma-1)/self.gamma))) - self.A_e/self.A_t

    def compute_magnitude(self, time, y):
        M_p, b = y

        # Get the current burning surface
        try:
            self.S = self.geometry_model.burning_S(b)
        # Return 0 if the burn thickness is too high
        except ValueError:
            return np.asarray([0, 0]), 0

        # Stop the thrust if there is no surface burning anymore
        if self.S == 0:
            self.m_dot, self.p_c, self.p_e, self.V_e, self.F_T, self.r = 0, 0, 0, 0, 0, 0
        else:
            # Compute propulsion characteristics at this time step
            self.m_dot = self.S * self.r * self.rho_p                               # [kg/s] mass flow
            self.p_c = (self.c_real * self.rho_p * self.a * self.S / self.A_t)**(1/(1-self.n)) # [Pa] combustion chamber pressure
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

        # Return the thrust
        return np.asarray([-self.m_dot, self.r]), self.F_T

    def rk_4(self, dt=1e-3):
        t = 0 if self.last_t is None else self.last_t
        y = [self.M_p, self.b]
        k1, F_t_1 = self.compute_magnitude(t, y)
        k2, F_t_2 = self.compute_magnitude(t + dt/2, y + dt/2 * k1)
        k3, F_t_3 = self.compute_magnitude(t + dt/2, y + dt/2 * k2)
        k4, F_t_4 = self.compute_magnitude(t + dt, y + dt * k3)
        self.F_t = (F_t_1 + 2*F_t_2 + 2*F_t_3 + F_t_4) / 6
        self.last_t = t + dt
        y_der = (k1 + 2*k2 + 2*k3 + k4) / 6
        y_next = y + dt * y_der
        self.M_p, self.b = y_next
        return y_next, y_der

    def euler(self, dt=1e-3):
        t = 0 if self.last_t is None else self.last_t
        y = [self.M_p, self.b]
        k, F_t = self.compute_magnitude(t, y)
        self.F_t = F_t
        self.last_t = t + dt
        y_next = y + dt * k
        self.M_p, self.b = y_next
        return y_next, k

    def simulate_full_burn(self, dt=1e-3, use_rk4=True, *args, **kwargs):
        # x     = [M_p,   b]
        # x_dot = [m_dot, r]
        b_s, p_c_s, M_p_s = [], [], []
        y = [self.M_p, self.b]
        i = 0
        # Keep computing thrust until the magnitude settles to 0
        while self.last_t is None or np.sum(self.saved_magnitudes[-2:]) != 0:
            if use_rk4:
                y, der = self.rk_4(dt)
            else:
                y, der = self.euler(dt)
            if i % 500 == 0:
                self.saved_magnitudes.append(self.F_t)
                self.saved_burn_times.append(self.last_t)
                self.saved_m_dot_s.append(der[0])
                b_s.append(y[1])
                M_p_s.append(y[0])
            i += 1
        print(i)

        self.saved_burn_times.insert(0, 0), self.saved_magnitudes.insert(0, 0), self.saved_m_dot_s.insert(0, 0)
        
        self.magnitude_interpolator = interp1d(self.saved_burn_times, self.saved_magnitudes)
        self.m_dot_interpolator = interp1d(self.saved_burn_times, self.saved_m_dot_s)

        return self.saved_burn_times, self.saved_magnitudes, b_s, p_c_s, M_p_s

if __name__ == "__main__":
    from thrust.models.multi_fin import multi_fin_SRM
    SRM = multi_fin_SRM(R_o=0.24, R_i=0.175, N_f=20, w_f=0.02, L_f=0.05, L=1.05)
    SRM_thrust_model = SRM_thrust_rk4(SRM, A_t=0.065, epsilon=45)
    times, magnitudes, b_s, _, M_p_s = SRM_thrust_model.simulate_full_burn()
    print(magnitudes[-7:-2])
    print(b_s[-7:-2])
    print(M_p_s[-7:-2])