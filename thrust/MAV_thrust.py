import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

# Standard imports
import numpy as np
from scipy.interpolate import interp1d

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import environment

# Custom imports
from thrust.solid_thrust import SRM_thrust

class MAV_thrust:

    def __init__(self, ascent_model, angle, thrust_model, body_directions_y=0, body_directions_z=0, print_status=False):
        self.angle = angle
        self.thrust_model = thrust_model

        if type(self.thrust_model) == list and type(self.thrust_model[0]) in [float, int]:
            self.thrust_type = "constant"
            self.magnitude = self.thrust_model[0]
            self.Isp = self.thrust_model[1]
            self.burn_time = self.thrust_model[2]
        elif type(self.thrust_model) == SRM_thrust:
            self.thrust_type = "from_geometry"
            if print_status:
                print("Pre-simulating SRM burn...", end="\r")
            self.thrust_model.simulate_full_burn()
            if print_status:
                print("Pre-simulating SRM burn finished.")
            self.magnitude_function = self.thrust_model.magnitude_interpolator
            self.m_dot_function = self.thrust_model.m_dot_interpolator
            self.burn_time = self.thrust_model.saved_burn_times[-1]
        else:
            raise NotImplementedError("The thrust model `%s` does not correspond to anything implemented" % type(self.thrust_model))
        
        if type(body_directions_y) == list:
            body_direction_times = np.linspace(0, self.burn_time, len(body_directions_y))
            self.body_direction_y_function = interp1d(body_direction_times, body_directions_y)
        else:
            self.body_direction_y_function = lambda t: body_directions_y
        if type(body_directions_z) == list:
            body_direction_times = np.linspace(0, self.burn_time, len(body_directions_z))
            self.body_direction_z_function = interp1d(body_direction_times, body_directions_z)
        else:
            self.body_direction_z_function = lambda t: body_directions_z

        self.ascent_model = ascent_model
        self.coasting = False
        self.first_t = None
        self.last_t = None

    def compute_time_elapsed(self, time):
        if self.first_t is None:
            self.first_t = time
            self.time_elapsed = 0
        else:
            self.time_elapsed = time - self.first_t
        self.last_t = time

    def is_thrust_on(self, time):
        if self.coasting:
            return False # Thrust is off if vehicle is coasting

        if self.last_t != time:
            self.compute_time_elapsed(time)

        # Do not compute thrust anymore if specified burn time has elapsed
        if self.time_elapsed > self.burn_time:
            self.coasting = True # If thrust is off, remember that we enter a coasting phase
            return False
        return True

    def get_thrust_magnitude(self, time):
        if self.thrust_type == "constant":
            return self.magnitude
        elif self.thrust_type == "from_geometry":
            return self.magnitude_function(self.time_elapsed)

    def get_specific_impulse(self, time):
        if self.thrust_type == "constant":
            return self.Isp
        elif self.thrust_type == "from_geometry":
            return 0

    def get_mass_flow(self, time):
        if self.last_t != time:
            self.compute_time_elapsed(time)

        if self.thrust_type == "from_geometry":
            if self.time_elapsed > self.burn_time:
                return 0
            return -np.fabs(self.m_dot_function(self.time_elapsed))
        else:
            raise NotImplementedError("No mass flow model has been implemented in this case.")

    def get_body_fixed_thrust_direction(self):
        if self.time_elapsed > self.burn_time:
            return np.zeros((3,3))
        else:
            angle_y = self.body_direction_y_function(self.time_elapsed)
            angle_z = self.body_direction_z_function(self.time_elapsed)
            T_y = np.array([
                [np.cos(angle_y),   0,      -np.sin(angle_y)],
                [0,                 1,      0],
                [np.sin(angle_y),   0,      np.cos(angle_y)]
            ])
            T_z = np.array([
                [np.cos(angle_z),   np.sin(angle_z),    0],
                [-np.sin(angle_z),  np.cos(angle_z),    0],
                [0,                 0,                  1]
            ])
            body_fixed_thrust_direction = np.dot(T_y, T_z)

            # print(self.time_elapsed)
            # print(angle_y)
            # print(T_y)
            # print(angle_z)
            # print(T_z)
            # print(np.linalg.norm(body_fixed_thrust_direction, axis=0))
            # print(body_fixed_thrust_direction), input()
            
            return body_fixed_thrust_direction

    def get_thrust_orientation(self, time):
        # Get aerodynamic angle calculator
        aerodynamic_angle_calculator = self.ascent_model.current_body.flight_conditions.aerodynamic_angle_calculator
        ### Force thrust parallel to angle of attack
        #self.angle = aerodynamic_angle_calculator.get_angle(environment.angle_of_attack)
        ###
        # Set thrust in vertical frame and transpose it
        thrust_direction_vertical_frame = np.array([[0, np.sin(self.angle), - np.cos(self.angle)]]).T
        thrust_direction_vertical_frame = thrust_direction_vertical_frame / np.linalg.norm(thrust_direction_vertical_frame)
        # Retrieve rotation matrix from vertical to inertial frame from the aerodynamic angle calculator
        vertical_to_inertial_frame = aerodynamic_angle_calculator.get_rotation_matrix_between_frames(
            environment.AerodynamicsReferenceFrames.vertical_frame,
            environment.AerodynamicsReferenceFrames.inertial_frame)
        # Compute the thrust in the inertial frame
        thrust_inertial_frame = np.dot(vertical_to_inertial_frame,
                                    thrust_direction_vertical_frame)
        # Add contribution from the body fixed direction
        body_fixed_direction = self.get_body_fixed_thrust_direction()
        thrust_orientation = np.dot(body_fixed_direction, thrust_inertial_frame)
        return thrust_orientation