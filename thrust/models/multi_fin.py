import numpy as np

class multi_fin_SRM:

    def __init__(self, R_o, R_i, N_f, w_f, L_f, L, run_checks=True):
        self.run_checks = run_checks
        # Check geometry validity
        if self.run_checks:
            c_1 = w_f/2 < R_o - R_i
            c_2 = L_f < R_i
            c_3 = N_f * w_f < 2*np.pi*(R_i-L_f)
            if not c_1:
                raise ValueError("The half fin width 'w_f/2' has to be smaller than the distance between the inner and outer radius 'R_o - R_i'.")
            if not c_2:
                raise ValueError("The fin length 'L_f' must be smaller than the inner radius 'R_i'.")
            if not c_3:
                raise ValueError("The fins geometrically interfere with each other, since the sum of their width is larger that the innermost radius.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i
        self.N_f = N_f
        self.w_f = w_f
        self.L_f = L_f
        self.L = L

        self.P_fin = 2 * self.N_f * self.L_f

    def check_b(self, b):
        if self.run_checks:
            # Check the validity of the given burnt thickness
            if b >= (self.R_o - self.R_i):
                raise ValueError("The burnt thickness 'b' is too high.")

    def get_V_p(self):
        # Return the propellant volume
        return (np.pi*(self.R_o**2 - self.R_i**2) + self.N_f*self.w_f*self.L_f) * self.L

    def burning_S(self, b):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        # No burning area if the outer radius is exceeded
        if self.R_i + b >= self.R_o:
            return 0
        # Compute the burning perimeter of the tubular part
        P_tube = 2*np.pi * (self.R_i + b)
        # Compute the burning perimeter of the fins
        P_fin = self.P_fin if b < self.w_f/2 else 0
        # Compute the burning surface
        return (P_tube + P_fin) * self.L

    def plot_geometry(self, b=0, final_b=None, ax_in=None, add_title=True):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        # Take burnt thickness into account in the geometry
        R_o = self.R_o
        R_i = min(self.R_i + b, self.R_o)
        N_f = self.N_f
        w_f = max(self.w_f - 2*b, 0)
        L_f = max(self.L_f - b, 0)

        # Convert certesian coordinates to polar coordinates
        def make_polar(x, y, angle=0):
            theta = np.arctan(y/x) + angle
            r = np.sqrt(x**2+y**2)
            return theta, r

        from matplotlib import pyplot as plt
        thetas = np.linspace(0, 2*np.pi, 100)
        if ax_in is None:
            ax = plt.subplot(111, polar=True)
        else:
            ax = ax_in
        # Plot the outer ring
        ax.plot(thetas, np.ones(100)*self.R_o, color="black")
        # Compute the angle encompassing the half width of the fin
        thickness_rad = np.arctan(w_f/2/R_i)
        # Fill everything in red (as burning)
        ax.fill_between(thetas, 0, R_o, color="#f6902f")
        # Plot the fins if they have not burn out yet
        if b < self.w_f/2:
            for i in range(self.N_f):
                # Plot the side of the fins
                side_theta_a, side_r = make_polar(np.linspace(R_i-L_f, R_i, 100), np.ones(100)*(w_f/2), 2*np.pi/N_f*i)
                side_theta_b, side_r = make_polar(np.linspace(R_i-L_f, R_i, 100), -np.ones(100)*(w_f/2), 2*np.pi/N_f*i)
                # Compute correction due to curvature
                dr = np.fabs(R_i**2*(max(side_theta_a)-2*np.pi/N_f*i)/4**2-w_f**2/4)
                ax.plot(side_theta_a, side_r-dr, color="black")
                ax.plot(side_theta_b, side_r-dr, color="black")
                # Plot the top of the fins
                top_theta, top_r = make_polar(np.ones(100)*(R_i-L_f), np.linspace(-(w_f/2), (w_f/2), 100), 2*np.pi/N_f*i)
                ax.plot(top_theta, top_r-dr, color="black")
                # Plot the inner ring, leaving holes where fins are
                left_a, right_a = thickness_rad+2*np.pi/N_f*i, -thickness_rad+2*np.pi/N_f*(i+1)
                etas_a = np.linspace(left_a, right_a, 100)
                ax.plot(etas_a, np.ones(100)*(R_i), color="black")
                # Fill the grain area inside the fins
                ax.fill_betweenx(np.linspace(R_i-L_f+dr, R_i, 100), side_theta_a, side_theta_b, color="#bbb")
        else:
            # Otherwise, simply plot the tubular geometry
            ax.plot(thetas, np.ones(100)*(R_i), color="black")
        # Fill the grain area between the outer ring and the fins
        ax.fill_between(thetas, R_i, R_o, color="#bbb")
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i-self.L_f, self.R_i, self.R_o])
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
            type="multi fin",
            R_o=self.R_o,
            R_i=self.R_i,
            N_f=self.N_f,
            w_f=self.w_f,
            L_f=self.L_f,
            L=self.L
        )

    def __str__(self):
        return """Multi-fin SRM geometry with
        $L = %.3f$ [m], $R_o = %.3f$ [m], $R_i = %.3f$ [m], $N_f = %i$ [-], $w_f = %.3f$ [m], $L_f = %.3f$ [m]""" \
        % (self.L, self.R_o, self.R_i, self.N_f, self.w_f, self.L_f)