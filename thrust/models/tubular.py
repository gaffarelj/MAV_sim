import numpy as np

class tubular_SRM:

    def __init__(self, R_o, R_i, L):
        # Check geometry validity
        c_1 = R_i < R_o
        if not c_1:
            raise ValueError("The outer radius 'R_o' must be larger than the inner radius 'R_i'.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i
        self.L = L

    def check_b(self, b):
        # Check the validity of the given burnt thickness
        if b >= (self.R_o - self.R_i):
            raise ValueError("The burnt thickness 'b' is too high.")

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        # Compute the burning perimeter
        burning_P = 2 * np.pi * (self.R_i + b)
        # Compute the burning surface
        return burning_P * self.L

    def plot_geometry(self, b=0, save=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        from matplotlib import pyplot as plt
        thetas = np.arange(0, 2*np.pi, 0.01)
        ax = plt.subplot(111, polar=True)
        # Plot the inner ring
        ax.plot(thetas, [self.R_i+b]*len(thetas), color="black")
        # Plot the outer ring
        ax.plot(thetas, [self.R_o]*len(thetas), color="black")
        # Fill the grain area between the rings
        ax.fill_between(thetas, self.R_i+b, self.R_o, color="grey", alpha=0.5)
        ax.set_ylim([0,self.R_o*1.5])
        ax.set_yticks([self.R_i, self.R_o])
        plt.title("Tubular SRM geometry with:\n$R_i=%.2f$ [m], $R_o=%.2f$ [m], $L=%.2f$ [m], $b=%.2f$ [m]" % (self.R_i, self.R_o, self.L, b))
        plt.tight_layout()
        # Save or show the plot
        if save is not None:
            plt.savefig(save)
            plt.close()
        else:
            plt.show()