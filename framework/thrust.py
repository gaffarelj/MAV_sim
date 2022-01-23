import sys
# Add tudatpy path
sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

import numpy as np

from tudatpy.kernel.numerical_simulation import environment

class MAV_thrust:

    def __init__(self, ascent_model, angle, magnitude, Isp):
        self.angle = angle
        self.magnitude = magnitude
        self.Isp = Isp
        self.ascent_model = ascent_model
        self.coasting = False

    def is_thrust_on(self, time):
        if self.coasting:
            return False
        current_stage_dry_mass = self.ascent_model.stage_1_dry_mass if self.ascent_model.current_stage == 1 else self.ascent_model.stage_2_dry_mass

        if self.ascent_model.current_body.mass > current_stage_dry_mass:
            return True
        else:
            self.coasting = True

    def get_thrust_magnitude(self, time):
        return self.magnitude

    def get_specific_impulse(self, time):
        return self.Isp

    def get_thrust_orientation(self, time):
        # Set thrust in vertical frame and transpose it
        thrust_direction_vertical_frame = np.array([[0, np.sin(self.angle), - np.cos(self.angle)]]).T
        # Get aerodynamic angle calculator
        aerodynamic_angle_calculator = self.ascent_model.current_body.flight_conditions.aerodynamic_angle_calculator
        # Retrieve rotation matrix from vertical to inertial frame from the aerodynamic angle calculator
        vertical_to_inertial_frame = aerodynamic_angle_calculator.get_rotation_matrix_between_frames(
            environment.AerodynamicsReferenceFrames.vertical_frame,
            environment.AerodynamicsReferenceFrames.inertial_frame)
        # Compute the thrust in the inertial frame
        thrust_inertial_frame = np.dot(vertical_to_inertial_frame,
                                    thrust_direction_vertical_frame)

        return thrust_inertial_frame