import numpy as np

class rod_and_tube_SRM:

    def __init__(self, R_o, R_mid, R_i, L, run_checks=True):
        self.run_checks = run_checks
        # Check geometry validity
        if self.run_checks:
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
        if self.run_checks:
            # Check the validity of the given burnt thickness
            if b > (self.R_o - self.R_mid) and b > self.R_i:
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

    def plot_geometry(self, b=0, final_b=None, ax_in=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        # Take burnt thickness into account in the geometry
        R_o = self.R_o
        R_mid = min(self.R_mid + b, self.R_o)
        R_i = max(self.R_i - b, 0)

        from matplotlib import pyplot as plt
        thetas = np.arange(0, 2*np.pi, 0.01)
        if ax_in is None:
            ax = plt.subplot(111, polar=True)
        else:
            ax = ax_in
        # Plot the inner ring
        ax.plot(thetas, [R_i]*len(thetas), color="black")
        # Plot the intermediate ring
        ax.plot(thetas, [R_mid]*len(thetas), color="black")
        # Plot the outer ring
        ax.plot(thetas, [R_o]*len(thetas), color="black")
        # Fill everything in red (as burning)
        ax.fill_between(thetas, 0, R_o, color="#f6902f")
        # Fill the grain area inside the inner ring
        ax.fill_between(thetas, 0, R_i, color="#bbb", alpha=1)
        # Fill the grain area between the intermediate and outer rings
        ax.fill_between(thetas, R_mid, R_o, color="#bbb", alpha=1)
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_mid, self.R_o])
        # Add title
        ax.set_title("$b = %.3f$ [m]" % b)
        # Use a tight layout
        plt.tight_layout()
        # Return the axis used
        return ax

    def __str__(self):
        return """Rod and tube SRM geometry with
        $L = %.3f$ [m], $R_o = %.3f$ [m], $R_mid = %.3f$ [m], $R_i = %i$ [m]""" \
        % (self.L, self.R_o, self.R_mid, self.R_i)