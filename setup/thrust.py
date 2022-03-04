import sys
# Add tudatpy path
sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

import numpy as np

from tudatpy.kernel.numerical_simulation import environment

from thrust.solid_thrust import SRM_thrust

class MAV_thrust:

    def __init__(self, ascent_model, angle, thrust_model):
        self.angle = angle
        self.thrust_model = thrust_model
        if type(self.thrust_model) == list and type(self.thrust_model[0]) in [float, int]:
            self.thrust_type = "constant"
            self.magnitude = self.thrust_model[0]
            self.Isp = self.thrust_model[1]
            self.burn_time = self.thrust_model[2]
        elif type(self.thrust_model) == SRM_thrust:
            self.thrust_type = "from_geometry"
            print("Pre-running SRM burn simulation...", end="\r")
            self.thrust_model.simulate_full_burn()
            print("Pre-running SRM burn simulation finished.")
            self.magnitude_function = self.thrust_model.magnitude_interpolator
            self.m_dot_function = self.thrust_model.m_dot_interpolator
            self.burn_time = self.thrust_model.saved_burn_times[-1]
        else:
            raise NotImplementedError("The thrust model `%s` does not correspond to anything implemented" % type(self.thrust_model))
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

    def get_thrust_orientation(self, time):
        # Get aerodynamic angle calculator
        aerodynamic_angle_calculator = self.ascent_model.current_body.flight_conditions.aerodynamic_angle_calculator
        ### Force thrust parallel to angle of attack
        #self.angle = aerodynamic_angle_calculator.get_angle(environment.angle_of_attack)
        ###
        # Set thrust in vertical frame and transpose it
        thrust_direction_vertical_frame = np.array([[0, np.sin(self.angle), - np.cos(self.angle)]]).T
        # Retrieve rotation matrix from vertical to inertial frame from the aerodynamic angle calculator
        vertical_to_inertial_frame = aerodynamic_angle_calculator.get_rotation_matrix_between_frames(
            environment.AerodynamicsReferenceFrames.vertical_frame,
            environment.AerodynamicsReferenceFrames.inertial_frame)
        # Compute the thrust in the inertial frame
        thrust_inertial_frame = np.dot(vertical_to_inertial_frame,
                                    thrust_direction_vertical_frame)
        return thrust_inertial_frame