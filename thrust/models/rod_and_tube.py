import numpy as np

class rod_and_tube_SRM:

    def __init__(self, R_o, R_mid, R_i, L):
        # Check geometry validity
        c_1 = R_i < R_mid
        c_2 = R_mid < R_o
        if not c_1:
            raise ValueError("The intermediate radius 'R_mid' must be larger than the inner radius 'R_i'.")
        if not c_2:
            raise ValueError("The outer radius 'R_o' must be larger than the intermediate radius 'R_mid'.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_mid = R_mid
        self.R_i = R_i
        self.L = L

    def check_b(self, b):
        # Check the validity of the given burnt thickness
        if b >= (self.R_o - self.R_mid) and b >= self.R_i:
            raise ValueError("The burnt thickness 'b' is too high.")

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        P_inner, P_outer = 0, 0
        # Compute the outer burning perimeter
        if self.R_mid + b < self.R_o:
            P_outer = 2 * np.pi * (self.R_mid + b)
        # Compute the inner burning perimeter
        if self.R_i - b > 0:
            P_inner = 2 * np.pi * (self.R_i - b)
        # Compute the burning area
        return (P_inner + P_outer) * self.L

    def plot_geometry(self, b=0, save=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        from matplotlib import pyplot as plt
        thetas = np.arange(0, 2*np.pi, 0.01)
        ax = plt.subplot(111, polar=True)
        # Plot the inner ring (only if it has not burned out yet)
        if self.R_i-b > 0:
            ax.plot(thetas, [self.R_i-b]*len(thetas), color="black")
        # Plot the intermediate ring (only if it has not burned out yet)
        if self.R_mid+b < self.R_o:
            ax.plot(thetas, [self.R_mid+b]*len(thetas), color="black")
        # Plot the outer ring (only if it has not burned out yet)
        if self.R_mid+b < self.R_o:
            ax.plot(thetas, [self.R_o]*len(thetas), color="black")
        # Fill everything in red (as burning)
            ax.fill_between(thetas, 0, self.R_o, color="#f6902f")
        # Fill the grain area inside the inner ring (only if it has not burned out yet)
        if self.R_i-b > 0:
            ax.fill_between(thetas, 0, self.R_i-b, color="#bbb", alpha=1)
        # Fill the grain area between the intermediate and outer rings (only if it has not burned out yet)
        if self.R_mid+b < self.R_o:
            ax.fill_between(thetas, self.R_mid+b, self.R_o, color="#bbb", alpha=1)
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_mid, self.R_o])
        plt.title(
            """Rod and tube SRM geometry with:
            $R_i=%.2f$ [m], $R_{mid}=%.2f$ [m], $R_o=%.2f$ [m],
            $L=%.2f$ [m], $b=%.2f$ [m]""" % \
             (self.R_i, self.R_mid, self.R_o, self.L, b))
        plt.tight_layout()
        # Save or show the plot
        if save is not None:
            plt.savefig(save)
            plt.close()
        else:
            plt.show()