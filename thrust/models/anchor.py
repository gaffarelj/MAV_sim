import numpy as np

class anchor_SRM:

    def __init__(self, R_o, R_i, N_a, w, r_f, delta_s, L):
        # Check geometry validity
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
        # Check the validity of the given burnt thickness
        if b >= max(self.w, self.w_t, self.w_ds):
            raise ValueError("The burnt thickness 'b' is too high.")

    def burning_S(self, b):
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

        return (P1+P2+P3+P4+P5+P6+P7) * 2 * self.N_a * self.L

    def plot_geometry(self, b=0, save=None, ax_in=None):
        # First, check the validity of the given burnt thickness
        self.check_b(b)

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
        ax.plot(thetas, [self.R_o]*len(thetas), color="black")
        # Compute the angle encompassing the half width of the anchor shank
        thickness_rad = np.arctan(self.delta_s/2/self.R_i)
        # Compute the angle encompassing the half width of the anchor spoke
        anchor_angle = np.pi/self.N_a - np.arctan((self.r_f+self.w)/(self.R_i+2*self.w+self.r_f))
        # Fill everything in grey (as propellant)
        ax.fill_between(thetas, 0, self.R_o, color="#bbb")
        for i in range(self.N_a):
            # Plot the side of the anchor
            side_theta_a, side_r = make_polar(np.linspace(self.R_i, self.R_i+2*self.w, 100), np.ones(100)*(self.delta_s/2), 2*np.pi/self.N_a*i)
            side_theta_b, side_r = make_polar(np.linspace(self.R_i, self.R_i+2*self.w, 100), -np.ones(100)*(self.delta_s/2), 2*np.pi/self.N_a*i)
            # Compute correction due to curvature
            dr = np.fabs(self.R_i**2*(max(side_theta_a)-2*np.pi/self.N_a*i)/4**2-self.delta_s**2/4)
            ax.plot(side_theta_a, side_r-dr, color="black")
            ax.plot(side_theta_b, side_r-dr, color="black")
            # Fill the burning area inside the fins
            ax.fill_betweenx(np.linspace(self.R_i+dr, self.R_i+2*self.w+b, 100), side_theta_a, side_theta_b, color="#f6902f")
            # Plot the inner circle with holes where the anchor is
            left_a, right_a = thickness_rad+2*np.pi/self.N_a*i, -thickness_rad+2*np.pi/self.N_a*(i+1)
            etas_a = np.linspace(left_a, right_a, 100)
            ax.plot(etas_a, np.ones(100)*(self.R_i), color="black")
            # Fill the burning area in the inner circle
            ax.fill_between(thetas, 0, self.R_i, color="#f6902f")
            # Plot the innermost side of the anchor
            left_b1, right_b1 = side_theta_a[-1], anchor_angle+2*np.pi/self.N_a*i
            etas_b1 = np.linspace(left_b1, right_b1, 100)
            ax.plot(etas_b1, np.ones(100)*(self.R_i+2*self.w), color="black")
            left_b2, right_b2 = side_theta_b[-1], -anchor_angle+2*np.pi/self.N_a*i
            etas_b2 = np.linspace(left_b2, right_b2, 100)
            ax.plot(etas_b2, np.ones(100)*(self.R_i+2*self.w), color="black")
            # Plot the outermost side of the anchor
            etas_b3 = np.linspace(right_b1, right_b2, 100)
            ax.plot(etas_b3, np.ones(100)*(self.R_o-self.w), color="black")
            # Fill the burning area in the anchor
            ax.fill_between(etas_b3, self.R_i+2*self.w, self.R_o-self.w, color="#f6902f")
            # Plot the anchor fillets
            fillet_center_x_a = self.R_i+2*self.w+self.r_f
            fillet_center_x_b = self.R_o-self.w-self.r_f
            circle_theta_a, circle_r_a = make_polar(fillet_center_x_a + np.cos(np.linspace(np.pi/2, np.pi, 100)) * self.r_f, \
                np.sin(np.linspace(np.pi/2, np.pi, 100)) * self.r_f, anchor_angle + 2*np.pi/self.N_a*i)
            ax.plot(circle_theta_a, circle_r_a, color="black")
            circle_theta_b, circle_r_b = make_polar(fillet_center_x_b + np.cos(np.linspace(0, np.pi/2, 100)) * self.r_f, \
                np.sin(np.linspace(0, np.pi/2, 100)) * self.r_f, anchor_angle + 2*np.pi/self.N_a*i)
            ax.plot(circle_theta_b, circle_r_b, color="black")
            circle_theta_c, circle_r_c = make_polar(fillet_center_x_a + np.cos(np.linspace(np.pi, 3*np.pi/2, 100)) * self.r_f, \
                np.sin(np.linspace(np.pi, 3*np.pi/2, 100)) * self.r_f, -anchor_angle + 2*np.pi/self.N_a*i)
            ax.plot(circle_theta_c, circle_r_c, color="black")
            circle_theta_d, circle_r_d = make_polar(fillet_center_x_b + np.cos(np.linspace(3*np.pi/2, 2*np.pi, 100)) * self.r_f, \
                np.sin(np.linspace(3*np.pi/2, 2*np.pi, 100)) * self.r_f, -anchor_angle + 2*np.pi/self.N_a*i)
            ax.plot(circle_theta_d, circle_r_d, color="black")
            # Plot fillet straight part
            fillet_line_a, fillet_line_b = [circle_theta_a[0], circle_theta_b[-1]], [circle_theta_c[-1], circle_theta_d[0]]
            ax.plot(fillet_line_a, [circle_r_a[0], circle_r_b[-1]], color="black")
            ax.plot(fillet_line_b, [circle_r_c[-1], circle_r_d[0]], color="black")
            # Fill the burning area from the fillet
            ax.fill_between(circle_theta_a, circle_r_a, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_b, circle_r_b, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_c, circle_r_c, np.ones(100)*fillet_center_x_a, color="#f6902f")
            ax.fill_between(circle_theta_d, circle_r_d, np.ones(100)*fillet_center_x_a, color="#f6902f")
        # Set plot limits and ticks
        ax.set_ylim([0,self.R_o*1.25])
        ax.set_yticks([self.R_i, self.R_i+2*self.w, self.R_o-self.w, self.R_o])
        # Use a tight layout
        plt.tight_layout()
        # Save or show the plot
        if save is not None:
            plt.savefig(save)
            plt.close()
        elif ax_in is None:
            plt.show()