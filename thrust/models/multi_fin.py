import numpy as np

class multi_fin_SRM:

    def __init__(self, R_o, R_i, N_f, w_f, L_f, L):
        # Check geometry validity
        c_1 = w_f/2 < R_o - R_i
        c_2 = L_f < R_i
        c_3 = np.sin(np.pi/N_f) * (w_f + 2 * (R_i - L_f)) >= 0
        if not c_1:
            raise ValueError("The half fin width 'w_f/2' has to be smaller than the distance between the inner and outer radius 'R_o - R_i'.")
        if not c_2:
            raise ValueError("The fin length 'L_f' must be smaller than the outer radius 'R_i'.")
        if not c_3:
            raise ValueError("The fin geometrically interfere with each other, since the following expression is negative: sin(pi/N_f) * (w_f + 2 * (R_i - L_f))")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i
        self.N_f = N_f
        self.w_f = w_f
        self.L_f = L_f
        self.L = L

    def check_b(self, b):
        # Check the validity of the given burnt thickness
        if b >= (self.R_o - self.R_i):
            raise ValueError("The burnt thickness 'b' is too high.")

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        # Compute the burning perimeter of the tubular part
        P_tube = 2*np.pi * (self.R_i + b)
        # Compute the burning perimeter of the fins
        P_fin = 0
        if b < self.w_f/2:
            P_fin = 2 * self.N_f * (self.L_f - b)
        # Compute the burning surface
        return (P_tube + P_fin) * self.L

    def plot_geometry(self, b=0, save=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        from matplotlib import pyplot as plt
        thetas = np.linspace(0, 2*np.pi, 100)
        ax = plt.subplot(111, polar=True)
        # Plot the outer ring
        ax.plot(thetas, [self.R_o]*len(thetas), color="black")
        # Compute the angle encompassing the half width of the fin
        thickness_rad = np.arcsin(self.w_f/2*1.25) - b
        # Plot the fins if they have not burn out yet
        if b < self.w_f/2:
            for i in range(self.N_f):
                # Plot the inner ring, leaving holes where fins are
                left_a, right_a = thickness_rad+2*np.pi/self.N_f*i, -thickness_rad+2*np.pi/self.N_f*(i+1)
                etas_a = np.linspace(left_a, right_a, 100)
                ax.plot(etas_a, [self.R_i+b]*len(etas_a), color="black")
                # Plot the side of the fins
                left_b, right_b = -thickness_rad+2*np.pi/self.N_f*i, thickness_rad+2*np.pi/self.N_f*i
                ax.vlines(left_b, self.R_i-self.L_f+b, self.R_i+b, color="black")
                ax.vlines(right_b, self.R_i-self.L_f+b, self.R_i+b, color="black")
                # Fill the grain area inside the fins
                ax.fill_betweenx(np.linspace(self.R_i-self.L_f+b, self.R_i+b, 100), left_b, right_b, color="grey", alpha=0.5)
                # Plot the top of the fins
                etas_b = np.linspace(left_b, right_b, 100)
                ax.plot(etas_b, [self.R_i-self.L_f+b]*len(etas_b), color="black")
        else:
            # Otherwise, simply plot the tubular geometry
            ax.plot(thetas, [self.R_i+b]*len(thetas), color="black")
        # Fill the grain area between the outer ring and the fins
        ax.fill_between(thetas, self.R_i+b, self.R_o, color="grey", alpha=0.5)
        ax.set_ylim([0,self.R_o*1.5])
        ax.set_yticks([self.R_i-self.L_f, self.R_i, self.R_o])
        plt.title("Tubular SRM geometry with:\n$R_i=%.2f$ [m], $R_o=%.2f$ [m], $L=%.2f$ [m], $b=%.2f$ [m]" % (self.R_i, self.R_o, self.L, b))
        plt.tight_layout()
        # Save or show the plot
        if save is not None:
            plt.savefig(save)
            plt.close()
        else:
            plt.show()