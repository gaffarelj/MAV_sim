from ctypes import c_int16
import numpy as np
from matplotlib import pyplot as plt

# Define geometry
N_a = 3
R_o = 1
R_i = 0.25
w = 0.2
r_f = 0.035
delta_s = 0.12

# Check conditions
c_1 = 0 < R_i and R_i < R_o
c_2 = 0 < w and w < (R_o - R_i) / 3
c_3 = 0 < r_f and r_f < (R_o - 3 * w - R_i) / 2
c_4 = 2 <= N_a
c_5 = 0 < delta_s and delta_s < 2* R_i * np.sin(np.pi/N_a)
c_6 = np.arcsin( (delta_s + 2 * w)/(2 * (R_i + w)) ) + np.arcsin( (r_f + w)/(R_i + 2 * w + r_f) ) < np.pi/N_a
print(c_1, c_2, c_3, c_4, c_5, c_6)

# Compute burning perimeter
w_t = np.sqrt((w + r_f)**2 + 2 * R_o**2 - 2 * R_o * np.sqrt(R_o**2 - 2 * R_o * (w + r_f)) + w + r_f ) - r_f
w_ds = (R_i**2 + ( 2 * w + r_f - np.sqrt(R_i**2 + 2 * R_i * (2 * w + r_f) + 3 * w**2 + 2 * r_f * w) ) * R_i + 2 * w * (w + r_f)) / \
    (np.sqrt(R_i**2 + 2 * R_i * (2 * w + r_f) + 3 * w**2 + 2 * r_f * w) + r_f - R_i)

P_s = []
b_s = np.arange(0, max(w, w_t, w_ds), 0.0001)
for b in b_s:
    # Perimeter 1
    if 0 <= b and b <= w:
        P1 = (R_i + b) * (np.pi/N_a - np.arcsin((delta_s+2*b)/(2*(R_i+b))))
    elif w < b and b <= w_ds:
        P1 = (R_i+b) * ( np.arcsin( (r_f+w)/(R_i+2*w+r_f) ) - np.arccos( (R_i**2+(2*w+b+r_f)*R_i+(2*w-b)*r_f+2*w**2)/((R_i+2*w+r_f)*(R_i+b)) ) )
    else:
        P1 = 0

    # Perimeter 2
    if b <= w:
        P2 = np.sqrt((R_i+2*w-b)**2-(delta_s/2+b)**2) - np.sqrt((R_i+b)**2-(delta_s/2+b)**2)
    else:
        P2 = 0

    # Perimeter 3
    if b <= w:
        P3 = (R_i+2*w-b)*( np.pi/N_a - np.arcsin( (delta_s+2*b)/(2*(R_i+2*w-b)) ) - np.arcsin( (r_f+w)/(R_i+2*w+r_f) ) )
    else:
        P3 = 0

    # Perimeter 4
    if b <= w:
        P4 = (r_f+b) * np.arccos( (r_f+w)/(R_i+2*w+r_f) )
    elif b <= w_ds:
        P4 = (r_f+b) * (np.arccos( (r_f+w)/(R_i+2*w+r_f) ) -np.arccos( (r_f+w)/(r_f+b) ) - np.arccos( ((R_i+2*w+r_f)**2 + (r_f+b)**2 - (R_i+b)**2)/(2*(R_i+2*w+r_f)*(r_f+b)) ) )
    else:
        P4 = 0

    # Perimeter 5
    if b <= w:
        P5 = np.sqrt((R_o-w-r_f)**2 - (r_f+w)**2) - np.sqrt((R_i+2*w+r_f)**2 - (r_f+w)**2)
    else:
        P5 = 0

    # Perimeter 6
    if b <= w:
        P6 = (r_f+b)* (np.pi/2 + np.arcsin( (r_f+w)/(R_o-w-r_f) ))
    elif b <= w_t:
        P6 = (r_f+b)* (np.arccos( ((R_o-w-r_f)**2 + (r_f+b)**2 - R_o**2)/(2*(R_o-w-r_f)*(r_f+b)) ) - np.arccos((r_f+w)/(r_f+b)) - np.arccos((r_f+w)/(R_o-w-r_f)))
        P6 = max(0, P6)
    else:
        P6 = 0

    # Perimeter 7
    if b <= w:
        P7 = (R_o-w+b) * (np.pi/N_a - np.arcsin( (r_f+w)/(R_o-w-r_f) ))
    else:
        P7 = 0

    P_s.append([P1, P2, P3, P4, P5, P6, P7])

P_s = np.array(P_s)

# for i in range(7):
#     plt.plot(b_s, P_s[:,i], label="P%i" % (i+1))
plt.plot(b_s, np.sum(P_s, axis=1)*2*N_a, label="Total")

plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning perimeter [m]")
# plt.legend()
plt.grid(), plt.tight_layout()
plt.show()
    

if False:
    ## Define geometry (http://dx.doi.org/10.5937/str1802048T, https://arc-aiaa-org.tudelft.idm.oclc.org/doi/pdf/10.2514/6.2008-4697)
    # R1 < R2 < R3 < R4 [m]
    R1, R2, R3, R4 = 0.1, 0.4, 0.45, 0.8
    # R1 < D < R2 [m]
    D = 0.25
    # S < R1 [m]
    S = 0.05
    # Length [m]
    L = 3

    def plot_geometry(R1, R2, R3, R4, D, S):
        def get_xy(R, angle):
            return R * np.sin(angle), R * np.cos(angle)

        def do_append(x, y, X, Y, condition=True):
            if condition:
                X.append(x), Y.append(y)
            else:
                X.append(np.nan), Y.append(np.nan)

        X1, Y1, X2, Y2, X3, Y3, X4, Y4 = [], [], [], [], [], [], [], []
        s2_a = np.inf
        for theta in np.linspace(0, 2*np.np.pi, 300):
            x1, y1 = get_xy(R1, theta)
            do_append(x1, y1, X1, Y1, np.fabs(x1) >= S/2)
            x2, y2 = get_xy(R2, theta)
            if R1 < R2:
                do_append(x2, y2, X2, Y2, np.fabs(y2) >= D and np.fabs(x2) >= S/2)
                if np.fabs(y2) < D:
                    s2_a = min(np.fabs(x2), s2_a)
            else:
                if np.fabs(y2) < D:
                    s2_a = min(np.fabs(x1), s2_a)
            x3, y3 = get_xy(R3, theta)
            do_append(x3, y3, X3, Y3, np.fabs(y3) >= D)
            x4, y4 = get_xy(R4, theta)
            do_append(x4, y4, X4, Y4)
        plt.plot(X1, Y1), plt.plot(X2, Y2), plt.plot(X3, Y3), plt.plot(X4, Y4)
        s2_b = s2_a + R3 - R2
        for s1 in [-1, +1]:
            for s2 in [-1, +1]:
                if R2 > R1:
                    plt.plot([s1*S/2, s1*S/2], [s2*R1, s2*R2], color="black")
                plt.plot([s1*s2_a, s1*s2_b], [s2*D, s2*D], color="black")
        plt.xlim([-1.0, 1.0]), plt.ylim([-1, 1])
        plt.axis('equal'), plt.grid()

    t = 0
    A_s = []
    while R3 < R4:
        # Phase 1
        if t <= (R2-R1)/2:
            P1 = (R1+t) * np.arccos((S/2+t)/(R1+t))
            P2 = (R2-t) * np.sin( np.arccos((S/2+t)/(R2-t)) ) - (R1+t) * np.sin( np.arccos((S/2+t)/(R1+t)) )
            P3 = (R2-t) * (np.np.pi/2 - np.arcsin((D-t)/(R2-t)) - np.arcsin((S/2+t)/(R2-t)))
            P4 = (R3+t) * np.cos( np.arcsin((D-t)/R3+t) ) - (R2-t)*np.cos( np.arcsin((D-t)/(R2-t)) )
            P5 = (R3+t) * (np.np.pi/2 - np.arcsin((D-t)/(R3+t)))
            P = 4*(P1+P2+P3+P4+P5)
            A = P*L
            R1_i, R3_i = R1, R3
        # Phase 2
        elif t <= D - (R2-R1)/2:
            D_i = D - (R2-R1)/2
            P1 = (R1_i+t) * np.arcsin((D_i-t)/(R1_i+t))
            P2 = (R3_i+t) * np.cos( np.arcsin((D_i-t)/(R3_i+t)) ) - (R1_i+t) * np.cos( np.arcsin((D_i-t)/(R1_i+t)) )
            P3 = (R3_i+t) * (np.np.pi/2 - np.arcsin((D_i-t)/(R3_i+t)))
            P = 4*(P1+P2+P3)
            A = P*L
        # Phase 3
        elif t <= R4-R3-D:
            pass
        R1 += t
        R2 -= t
        R3 += t
        S += 2*t
        D -= t
        # plot_geometry(R1, R2, R3, R4, D, S)
        # plt.draw()
        # plt.pause(0.0005)
        # plt.clf()
        t += 0.0001
        A_s.append(A)
    plt.plot(np.linspace(0, t, len(A_s)), A_s)
    plt.xlabel("Burnt thickness [m]"), plt.ylabel("Burning area [m$^2$]")
    plt.grid()
    plt.tight_layout()
    plt.show()


if False:
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
    dS = 0.5    # [m2] difference of surface burned betweet t_i and t_i-1
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
