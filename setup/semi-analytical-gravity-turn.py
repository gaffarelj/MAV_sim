from matplotlib import pyplot as plt
import numpy as np
import sys


use_grav_turn = True
use_Drag = True
plot_rho_scale = False

phase = "2" # 1, coasting, 2

for phase in ["1", "coasting", "2"]:

        ### Initial parameters
        if phase == "1":
                V0 = 100
                h0 = 0
                M_P = 216
                M_0 = 385
                F_T = 9854
                c_eff = 2874.33
                delta = 30*np.pi/180
                if not use_grav_turn:
                        delta = 0
                gamma_0 = np.pi/2 - delta
        elif phase == "coasting":
                V0 = V_s[-1]
                h0 = Z_s[-1]
                M_P = 0
                M_0 = 385-216
                F_T = 0
                c_eff = 0
                gamma_0 = gamma
        elif phase == "2":
                V0 = V_s[-1]
                h0 = Z_s[-1]
                M_P = 54
                M_0 = 99
                F_T = 6937
                c_eff = 2766.42
                delta = np.pi/2
                gamma_0 = np.pi/2 - delta
        else:
                raise NotImplementedError

        g_0 = 9.80665
        g = 3.721
        psi_0 = F_T / (M_0*g_0)
        Lambda = M_0/(M_0-M_P)
        CD = 1.1
        S_ref = 0.144
        dt = 0.001

        def rho_expo(h):
                if h <= 25e3:
                        rho_0 = 0.0159
                        h_s = 11049
                else:
                        rho_0 = 0.0525
                        h_s = 7295
                return rho_0 * np.exp(-h/h_s)

        if plot_rho_scale:
                hs = np.arange(0, 120e3, 100)
                rhos= [rho_expo(h) for h in hs]
                plt.plot(hs/1e3, rhos)
                plt.yscale("log")
                plt.ylabel("Density [kg/m3]")
                plt.xlabel("Altitude [km]")
                plt.grid()
                plt.tight_layout()
                plt.show()

        ### Constant Thrust T

        # Compute burn time
        try:
                t_b_T = c_eff/g_0/psi_0 * (1-1/Lambda)
                print(phase, ": t_b_T =", t_b_T, "s")
                m_dt = psi_0 * g_0 / c_eff
                mdot = F_T / c_eff
        except ZeroDivisionError:
                t_b_T = 0
                m_dt = 0
                mdot = 0
        # Initialise the variables
        time, M_M0, gamma = 0, 1, gamma_0
        t0 = 0
        M = M_0
        X, Z = 0, h0
        V = V0
        # Compute the constant mass flow and Thrust
        T = c_eff * m_dt
        # Create empty arrays to store the results
        if phase == "1":
                time_s, V_s, M_M0_s, Z_s, X_s = [], [], [], [], []
                M_s = []
        else:
                time = time_s[-1]
                t0 = time
        V_loss = []
        V_loss_drag = []
        Vz = 0.01
        # Keep going until next step is above burn time
        print("Simulating constant Thrust T...")
        while (phase == "coasting" and Vz > 0) or (phase != "coasting" and time-t0+dt <= t_b_T):
                time += dt
                # Compute the new gamma over time
                dgamma_dt = -g*np.cos(gamma)/V
                # Compute the new mass and mass ratio
                M_M0 -= dt * m_dt
                M -= dt * mdot
                # Compute the new Velocity over time and velocity
                a = T/M_M0
                a_grav = g*np.sin(gamma)
                a_drag = 0
                if use_Drag:
                        rho = rho_expo(Z)
                        q = 0.5*rho*V**2
                        a_drag = CD * S_ref * q / M
                        M_s.append(M)
                V_loss_drag.append(a_drag*dt)
                dV_dt = a - a_grav - a_drag
                V_loss.append(a_grav*dt)
                V += dV_dt * dt
                # Compute the new velocities and positions
                Vx = V*np.cos(gamma)
                Vz = V*np.sin(gamma)
                X += Vx * dt
                Z += Vz * dt
                # Compute the new gamma
                gamma += dgamma_dt * dt
                # Store the new data
                V_s.append(V), M_M0_s.append(M_M0), X_s.append(X), Z_s.append(Z)
                time_s.append(time)

        grav_V_loss = sum(V_loss)
        drag_V_loss = sum(V_loss_drag)
        tot_dV = V_s[-1]-V0
        dV_thrust = tot_dV+grav_V_loss+drag_V_loss

        print("Phase is", phase)
        if phase == "coasting":
                print("Total coasting time:", time-t0, "s")
        print("Total dV of %.2f m/s (thrust of %.2f m/s, gravity loss of %.2f m/s, drag loss of %.2e m/s)" % (tot_dV, dV_thrust, grav_V_loss, drag_V_loss))
        print("Final altitude of %.5f km" % (Z_s[-1]/1e3))
        print("Final velocity of %.3f m/s" % V_s[-1])
        print("Final gamma of %.3f deg" % np.rad2deg(gamma))

### Plot the results
plt.figure(figsize=(9,3.5))
plt.subplot(121)
plt.plot(np.asarray(time_s), np.asarray(Z_s)/1e3)
plt.xlabel("Time [s]")
plt.ylabel("Altitude [km]")
plt.grid()
plt.subplot(122)
plt.plot(time_s, V_s)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.grid()
plt.tight_layout()
plt.savefig(sys.path[0]+"/../plots/setup/semi_analytical.pdf")
