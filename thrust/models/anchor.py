import numpy as np

class anchor_SRM:

    def __init__(self, R_o, R_i, N_a, w, r_f, delta_s, L, run_checks=True):
        self.run_checks = run_checks
        # Check geometry validity
        if self.run_checks:
            c_1 = 0 < R_i and R_i < R_o
            c_2 = 0 < w and w < (R_o - R_i) / 3
            c_3 = 0 < r_f and r_f < (R_o - 3 * w - R_i) / 2
            c_4 = 2 <= N_a
            c_5 = 0 < delta_s and delta_s < 2* R_i * np.sin(np.pi/N_a)
            c_6 = np.arcsin( (delta_s + 2 * w)/(2 * (R_i + w)) ) + np.arcsin( (r_f + w)/(R_i + 2 * w + r_f) ) < np.pi/N_a
            if not c_1:
                raise ValueError("The inner radius 'R_i' must be smaller than the outer radius 'R_o'.")
            if not c_2:
                raise ValueError("The anchor spacing 'w' must be at most a third of the distance between 'R_o' and 'R_i'.")
            if not c_3:
                raise ValueError("The fillet radius 'r_f' must be smaller than the remaining space when removing the inner radius 'R_i' and 3 times the anchor spacing 'w' from the outer radius 'R_o'.")
            if not c_4:
                raise ValueError("There should be at least two anchor spokes.")
            if not c_5:
                raise ValueError("The width of the anchor shank must allow for the inner arc defined by 'R_i' to be higher than a point.")
            if not c_6:
                raise ValueError("With this geometry, the anchors will be attached together while the centre of the SRM will still have propellant, which is not valid.")
        # Save the parameters of the geometry
        self.R_o = R_o
        self.R_i = R_i
        self.N_a = N_a
        self.w = w
        self.r_f = r_f
        self.delta_s = delta_s
        self.L = L
        # Compute the web thickness
        self.w_t = np.sqrt((w + r_f)**2 + 2 * R_o**2 - 2 * R_o * np.sqrt(R_o**2 - 2 * R_o * (w + r_f)) + w + r_f ) - r_f
        # Compute the web through the detached sliver
        self.w_ds = (R_i**2 + ( 2 * w + r_f - np.sqrt(R_i**2 + 2 * R_i * (2 * w + r_f) + 3 * w**2 + 2 * r_f * w) ) * R_i + 2 * w * (w + r_f)) / \
            (np.sqrt(R_i**2 + 2 * R_i * (2 * w + r_f) + 3 * w**2 + 2 * r_f * w) + r_f - R_i)

    def check_b(self, b):
        if self.run_checks:
            # Check the validity of the given burnt thickness
            if b >= max(self.w, self.w_t, self.w_ds):
                raise ValueError("The burnt thickness 'b' is too high.")

    def get_V_p(self):
        # Return the propellant volume
        A1 = 0.5*(np.pi/self.N_a - np.arcsin(self.delta_s/(2*self.R_i))) * self.R_i**2
        A2 = ((self.R_i+2*self.w)*self.R_i)/2*np.sin( np.arcsin(self.delta_s/(2*self.R_i)) - np.arcsin(self.delta_s/2 / (self.R_i+2*self.w)) )
        A3 = 0.5*np.arcsin(self.delta_s/(2*self.R_i+2*self.w))*(self.R_i+2*self.w)**2
        A4 = 0.5*(np.pi/self.N_a - np.arcsin((self.w+self.r_f)/(self.R_i+2*self.w+self.r_f))) * ((self.R_o-self.w)**2-(self.R_i+2*self.w)**2)
        A5 = 0.5*(np.arcsin((self.w+self.r_f)/(self.R_i+2*self.w+self.r_f)) - np.arcsin((self.w+self.r_f)/(self.R_o-self.w-self.r_f))) * (self.R_o-self.w)**2 - \
            0.5 * (self.R_o-self.w-self.r_f)*(self.R_i+2*self.w+self.r_f) * np.sin( np.arcsin((self.w+self.r_f)/(self.R_i+2*self.w+self.r_f)) - np.arcsin( (self.w+self.r_f)/(self.R_o-self.w-self.r_f)) )
        A6 = 0.5*np.arccos((self.w+self.r_f)/(self.R_i+2*self.w+self.r_f))*self.r_f**2+self.r_f*(np.sqrt((self.R_o-self.w-self.r_f)**2-(self.w+self.r_f)**2) - np.sqrt((self.R_i+2*self.w+self.r_f)**2-(self.w+self.r_f)**2))
        A7 = 0.5 * (np.pi/2 + np.arcsin((self.w+self.r_f)/(self.R_o-self.w-self.r_f)))*self.r_f**2
        return (np.pi * self.R_o**2 - 2 * self.N_a * (A1+A2+A3+A4+A5+A6+A7) ) * self.L

    def burning_S(self, b, return_perimeters=False):
        # First, check the validity of the given burnt thickness
        self.check_b(b)
        
        P1, P2, P3, P4, P5, P6, P7 = [0]*7
        # Perimeter 1
        if 0 <= b and b <= self.w:
            P1 = (self.R_i + b) * (np.pi/self.N_a - np.arcsin((self.delta_s+2*b)/(2*(self.R_i+b))))
        elif self.w < b and b <= self.w_ds:
            P1 = (self.R_i+b) * ( np.arcsin( (self.r_f+self.w)/(self.R_i+2*self.w+self.r_f) ) - \
                np.arccos( (self.R_i**2+(2*self.w+b+self.r_f)*self.R_i+(2*self.w-b)*self.r_f+2*self.w**2)/((self.R_i+2*self.w+self.r_f)*(self.R_i+b)) ) )
        # Perimeter 2
        if b <= self.w:
            P2 = np.sqrt((self.R_i+2*self.w-b)**2-(self.delta_s/2+b)**2) - np.sqrt((self.R_i+b)**2-(self.delta_s/2+b)**2)
        # Perimeter 3
        if b <= self.w:
            P3 = (self.R_i+2*self.w-b)*( np.pi/self.N_a - np.arcsin( (self.delta_s+2*b)/(2*(self.R_i+2*self.w-b)) ) -\
                np.arcsin( (self.r_f+self.w)/(self.R_i+2*self.w+self.r_f) ) )
        # Perimeter 4
        if b <= self.w:
            P4 = (self.r_f+b) * np.arccos( (self.r_f+self.w)/(self.R_i+2*self.w+self.r_f) )
        elif b <= self.w_ds:
            P4 = (self.r_f+b) * (np.arccos( (self.r_f+self.w)/(self.R_i+2*self.w+self.r_f) ) \
                -np.arccos( (self.r_f+self.w)/(self.r_f+b) ) - np.arccos( ((self.R_i+2*self.w+self.r_f)**2 + \
                (self.r_f+b)**2 - (self.R_i+b)**2)/(2*(self.R_i+2*self.w+self.r_f)*(self.r_f+b)) ) )
        # Perimeter 5
        if b <= self.w:
            P5 = np.sqrt((self.R_o-self.w-self.r_f)**2 - (self.r_f+self.w)**2) - np.sqrt((self.R_i+2*self.w+self.r_f)**2 - (self.r_f+self.w)**2)
        # Perimeter 6
        if b <= self.w:
            P6 = (self.r_f+b)* (np.pi/2 + np.arcsin( (self.r_f+self.w)/(self.R_o-self.w-self.r_f) ))
        elif b <= self.w_t:
            P6 = (self.r_f+b)* (np.arccos( ((self.R_o-self.w-self.r_f)**2 + (self.r_f+b)**2 - self.R_o**2)/(2*(self.R_o-self.w-self.r_f)*(self.r_f+b)) ) \
                - np.arccos((self.r_f+self.w)/(self.r_f+b)) - np.arccos((self.r_f+self.w)/(self.R_o-self.w-self.r_f)))
            #### TODO: check this (change w_t to w_ds ?)
            P6 = max(0, P6)
        # Perimeter 7
        if b <= self.w:
            P7 = (self.R_o-self.w+b) * (np.pi/self.N_a - np.arcsin( (self.r_f+self.w)/(self.R_o-self.w-self.r_f) ))

        A = (P1+P2+P3+P4+P5+P6+P7) * 2 * self.N_a * self.L
        if return_perimeters:
            return [P1, P2, P3, P4, P5, P6, P7], A
        return A

    def plot_geometry(self, b=0, final_b=None, ax_in=None, add_title=True):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

        # Take burnt thickness into account in the geometry
        R_o = self.R_o
        R_i = self.R_i + b
        N_a = self.N_a
        w = self.w - b
        r_f = self.r_f + b
        delta_s = self.delta_s + 2*b

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
        if final_b is None:
            _a = 1
        else:
            burn_progress = 0 if b == 0 else b/final_b
            _a = 1 if burn_progress < 0.2 else max((11-20*burn_progress)/7, 0)
        # Plot the outer ring
        ax.plot(thetas, [R_o]*len(thetas), color="black")
        # Compute the angle encompassing the half width of the anchor shank
        thickness_rad = np.arctan(delta_s/2/R_i)
        # Compute the angle encompassing the half width of the anchor spoke
        anchor_angle = np.pi/N_a - np.arctan((r_f+w)/(R_i+2*w+r_f))
        # Fill everything in grey (as propellant)
        ax.fill_between(thetas, 0, R_o, color="#bbb")
        for i in range(N_a):
            # Plot the side of the anchor
            side_theta_a, side_r = make_polar(np.linspace(R_i, R_i+2*w, 100), np.ones(100)*(delta_s/2), 2*np.pi/N_a*i)
            side_theta_b, side_r = make_polar(np.linspace(R_i, R_i+2*w, 100), -np.ones(100)*(delta_s/2), 2*np.pi/N_a*i)
            # Compute correction due to curvature
            dr = np.fabs(R_i**2*(max(side_theta_a)-2*np.pi/N_a*i)/4**2-delta_s**2/4)
            ax.plot(side_theta_a, side_r-dr, color="black", alpha=_a)
            ax.plot(side_theta_b, side_r-dr, color="black", alpha=_a)
            # Fill the burning area inside the fins
            ax.fill_betweenx(np.linspace(R_i-dr, R_i+2*w+r_f, 100), side_theta_a, side_theta_b, color="#f6902f")
            # Plot the inner circle with holes where the anchor is
            left_a, right_a = thickness_rad+2*np.pi/N_a*i, -thickness_rad+2*np.pi/N_a*(i+1)
            etas_a = np.linspace(left_a, right_a, 100)
            ax.plot(etas_a, np.ones(100)*(R_i), color="black", alpha=_a)
            # Fill the burning area in the inner circle
            ax.fill_between(thetas, 0, R_i, color="#f6902f")
            # Plot the innermost side of the anchor
            left_b1, right_b1 = side_theta_a[-1], anchor_angle+2*np.pi/N_a*i
            etas_b1 = np.linspace(left_b1, right_b1, 100)
            ax.plot(etas_b1, np.ones(100)*(R_i+2*w), color="black", alpha=_a)
            left_b2, right_b2 = side_theta_b[-1], -anchor_angle+2*np.pi/N_a*i
            etas_b2 = np.linspace(left_b2, right_b2, 100)
            ax.plot(etas_b2, np.ones(100)*(R_i+2*w), color="black", alpha=_a)
            # Plot the outermost side of the anchor
            etas_b3 = np.linspace(right_b1, right_b2, 100)
            ax.plot(etas_b3, np.ones(100)*(R_o-w), color="black", alpha=_a)
            # Fill the burning area in the anchor
            ax.fill_between(etas_b3, R_i+2*w, R_o-w, color="#f6902f")
            # Plot the anchor fillets
            fillet_center_x_a = R_i+2*w+r_f
            fillet_center_x_b = R_o-w-r_f
            circle_theta_a, circle_r_a = make_polar(fillet_center_x_a + np.cos(np.linspace(np.pi/2, np.pi, 100)) * r_f, \
                np.sin(np.linspace(np.pi/2, np.pi, 100)) * r_f, anchor_angle + 2*np.pi/N_a*i)
            ax.plot(circle_theta_a, circle_r_a, color="black", alpha=_a)
            circle_theta_b, circle_r_b = make_polar(fillet_center_x_b + np.cos(np.linspace(0, np.pi/2, 100)) * r_f, \
                np.sin(np.linspace(0, np.pi/2, 100)) * r_f, anchor_angle + 2*np.pi/N_a*i)
            ax.plot(circle_theta_b, circle_r_b, color="black", alpha=_a)
            circle_theta_c, circle_r_c = make_polar(fillet_center_x_a + np.cos(np.linspace(np.pi, 3*np.pi/2, 100)) * r_f, \
                np.sin(np.linspace(np.pi, 3*np.pi/2, 100)) * r_f, -anchor_angle + 2*np.pi/N_a*i)
            ax.plot(circle_theta_c, circle_r_c, color="black", alpha=_a)
            circle_theta_d, circle_r_d = make_polar(fillet_center_x_b + np.cos(np.linspace(3*np.pi/2, 2*np.pi, 100)) * r_f, \
                np.sin(np.linspace(3*np.pi/2, 2*np.pi, 100)) * r_f, -anchor_angle + 2*np.pi/N_a*i)
            ax.plot(circle_theta_d, circle_r_d, color="black", alpha=_a)
            # Plot fillet straight part
            fillet_line_a, fillet_line_b = [circle_theta_a[0], circle_theta_b[-1]], [circle_theta_c[-1], circle_theta_d[0]]
            ax.plot(fillet_line_a, [circle_r_a[0], circle_r_b[-1]], color="black", alpha=_a)
            ax.plot(fillet_line_b, [circle_r_c[-1], circle_r_d[0]], color="black", alpha=_a)
            # Fill the burning area from the fillet
            ax.fill_between(circle_theta_a, circle_r_a, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_b, circle_r_b, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_c, circle_r_c, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_d, circle_r_d, np.ones(100)*fillet_center_x_a, color="#f6902f")
        # Remove overflow
        ax.fill_between(thetas, self.R_o, 2*self.R_o, color="#fff")
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_i+2*self.w, self.R_o-self.w, self.R_o])
        # Add title
        if add_title:
            ax.set_title("$b = %.3f$ [m]" % b)
        # Use a tight layout
        plt.tight_layout()
        # Return the axis used
        return ax

    def __str__(self):
        return """Anchor SRM geometry with
        $L = %.3f$ [m], $R_o = %.3f$ [m], $R_i = %.3f$ [m], $N_a = %i$ [-], $w = %.3f$ [m], $r_f = %.3f$ [m], $\delta_s = %.3f$ [m]""" \
        % (self.L, self.R_o, self.R_i, self.N_a, self.w, self.r_f, self.delta_s)


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    # Define geometry from DOI: 10.2514/6.2008-4697
    anchor_verif = anchor_SRM(
        R_o=1.0,
        R_i=0.25,
        N_a=3,
        w=0.2,
        r_f=0.035,
        delta_s=0.12,
        L=1.0 # Set to 1 so the Area equals the Perimeter
    )

    # Plot the initial geometry
    anchor_verif.plot_geometry()
    plt.show()

    # Compute the burning area as a function of the burned distance
    P_s, b_s = [], []
    P_segments = []
    for b in np.arange(0, 1, 0.001):
        try:
            perimeters, A = anchor_verif.burning_S(b, True)
            P_segments.append(perimeters)
            P_s.append(A)
        except ValueError:
            break
        b_s.append(b)

    P_segments = np.array(P_segments).T

    # Plot the burning Perimeter as a function of the burned distance
    fig = plt.subplots(figsize=(8, 4))
    plt.plot(b_s, P_s)
    plt.xlim(0, 0.38), plt.ylim(0, 13)
    plt.xlabel("Burned distance [m]"), plt.ylabel("Burning Perimeter [m]")
    plt.grid()
    plt.tight_layout()
    plt.show()

    # Plot the burning perimeter segments as a function of the burned distance
    fig = plt.subplots(figsize=(8, 4))
    for i, perimeters in enumerate(P_segments):
        plt.plot(b_s, perimeters, label="Segment %i" % i)
    plt.xlabel("Burned distance [m]"), plt.ylabel("Burning Perimeter [m]")
    plt.xlim(0, 0.36), plt.ylim(0, 0.8)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
