import numpy as np

class spherical_SRM:

    def __init__(self, R_o, R_i, run_checks=True):
        self.run_checks = run_checks
        # Check geometry validity
        if self.run_checks:
            c_1 = R_i < R_o
            if not c_1:
                raise ValueError("The outer radius 'R_o' must be larger than the inner radius 'R_i'.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i

        # Assume exhaust tube of 50% inner radius
        self.R_e = 0.5*R_i

    def check_b(self, b):
        if self.run_checks:
            # Check the validity of the given burnt thickness
            if b >= (self.R_o - self.R_i):
                raise ValueError("The burnt thickness 'b' is too high.")

    def get_V_p(self):
        # Compute the exhaust tube volume
        V_e = np.pi * self.R_e ** 2 * (self.R_o - self.R_i)
        # Return the propellant volume
        return 4/3 * np.pi * (self.R_o**3 - self.R_i**3) - V_e

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        # No burning area if the outer radius is exceeded
        if self.R_i + b >= self.R_o:
            return 0
        # Compute the burning surface
        return 4 * np.pi * (self.R_i+b)**2 - np.pi * (self.R_e+b)**2 + 2*np.pi*(self.R_e+b)*(self.R_o - self.R_i-b)

    def plot_geometry(self, b=0, final_b=None, ax_in=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        # Take burnt thickness into account in the geometry
        R_o = self.R_o
        R_i = min(self.R_i + b, self.R_o)
        R_e = self.R_e + b

        # Convert certesian coordinates to polar coordinates
        def make_polar(x, y, angle=0):
            theta = np.arctan(y/x) + angle
            r = np.sqrt(x**2+y**2)
            return theta, r

        from matplotlib import pyplot as plt
        thetas = np.arange(0, 2*np.pi, 0.01)
        # Compute opening angle of exhaust tube at R_o and R_i
        open_R_o = np.arctan(R_e/R_o)
        open_R_i = np.arctan(R_e/R_i)*1.05
        if ax_in is None:
            ax = plt.subplot(111, polar=True)
        else:
            ax = ax_in
        # Plot the inner ring
        ax.plot(np.linspace(open_R_i, -open_R_i+2*np.pi, 100), np.ones(100)*R_i, color="black")
        # Plot the outer ring
        ax.plot(np.linspace(open_R_o, -open_R_o+2*np.pi, 100), np.ones(100)*R_o, color="black")
        # Plot the sides of the exhaust tube
        angle_exhaust = np.array([open_R_i, open_R_o])
        ax.plot(angle_exhaust, [R_i, R_o], color="black")
        ax.plot(-angle_exhaust, [R_i, R_o], color="black")
        # Fill everything in grey (as propellant)
        ax.fill_between(thetas, 0, R_o, color="#bbb")
        # Fill the burning area in the inner ring
        ax.fill_between(thetas, 0, R_i, color="#f6902f")
        # Fill the burning area in the exhaust tube
        ax.fill_betweenx(np.linspace(R_i, R_o, 2), -angle_exhaust, angle_exhaust, color="#f6902f")
        # Fill the remaining burning area outward from the exhaust tube
        ax.fill_between(np.linspace(-open_R_o, open_R_o, 100)*0.95, R_i, R_o, color="#f6902f")
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_o])
        # Add title
        ax.set_title("$b = %.3f$ [m]" % b)
        # Use a tight layout
        plt.tight_layout()
        # Return the axis used
        return ax

    def __str__(self):
        return """Spherical SRM geometry with
        $R_o = %.3f$ [m], $R_i = %i$ [m]""" \
        % (self.R_o, self.R_i)