import numpy as np

class tubular_SRM:

    def __init__(self, R_o, R_i, L, run_checks=True):
        self.run_checks = run_checks
        # Check geometry validity
        if self.run_checks:
            c_1 = R_i < R_o
            if not c_1:
                raise ValueError("The outer radius 'R_o' must be larger than the inner radius 'R_i'.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i
        self.L = L

    def check_b(self, b):
        if self.run_checks:
            # Check the validity of the given burnt thickness
            if b >= (self.R_o - self.R_i):
                raise ValueError("The burnt thickness 'b' is too high.")

    def get_V_p(self):
        # Return the propellant volume
        return np.pi*(self.R_o**2 - self.R_i**2) * self.L

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        # No burning area if the outer radius is exceeded
        if self.R_i + b >= self.R_o:
            return 0
        # Compute the burning perimeter
        burning_P = 2 * np.pi * (self.R_i + b)
        # Compute the burning surface
        return burning_P * self.L

    def plot_geometry(self, b=0, final_b=None, ax_in=None, add_title=True):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        # Take burnt thickness into account in the geometry
        R_o = self.R_o
        R_i = min(self.R_i + b, self.R_o)

        from matplotlib import pyplot as plt
        thetas = np.arange(0, 2*np.pi, 0.01)
        if ax_in is None:
            ax = plt.subplot(111, polar=True)
        else:
            ax = ax_in
        # Plot the inner ring
        ax.plot(thetas, [R_i]*len(thetas), color="black")
        # Plot the outer ring
        ax.plot(thetas, [R_o]*len(thetas), color="black")
        # Fill everything in red (as burning)
        ax.fill_between(thetas, 0, R_o, color="#f6902f")
        # Fill the grain area between the rings
        ax.fill_between(thetas, R_i, R_o, color="#bbb")
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_o])
        # Add title
        if add_title:
            ax.set_title("$b = %.3f$ [m]" % b)
        # Use a tight layout
        plt.tight_layout()
        # Return the axis used
        return ax

    def get_cpp_counterpart(self):
        from thrust.models.CPP.SRM_cpp import SRM_geometry
        return SRM_geometry(
            type="tubular",
            R_o=self.R_o,
            R_i=self.R_i,
            L=self.L
        )

    def __str__(self):
        return """Tubular SRM geometry with
        $L = %.3f$ [m], $R_o = %.3f$ [m], $R_i = %i$ [m]""" \
        % (self.L, self.R_o, self.R_i)