import sys
# Add tudatpy path
sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

import numpy as np

from tudatpy.kernel.numerical_simulation import environment

from thrust.solid_thrust import thrust as SRM_thrust

class MAV_thrust:

    def __init__(self, ascent_model, angle, thrust_model):
        self.angle = angle
        if type(thrust_model) == list and type(thrust_model[0]) in [float, int]:
            self.thrust_type = "constant"
            self.magnitude = thrust_model[0]
            self.Isp = thrust_model[1]
        elif type(thrust_model) == SRM_thrust:
            self.thrust_type = "from_geometry"
            self.magnitude_function = thrust_model.compute_magnitude
            self.Isp_function = thrust_model.get_Isp
        else:
            raise NotImplementedError("The thrust model `%s` does not correspond to anything implemented" % type(thrust_model))
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
        if self.thrust_type == "constant":
            return self.magnitude
        elif self.thrust_type == "from_geometry":
            return self.magnitude_function(time)

    def get_specific_impulse(self, time):
        if self.thrust_type == "constant":
            return self.Isp
        elif self.thrust_type == "from_geometry":
            return self.Isp_function(time)

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