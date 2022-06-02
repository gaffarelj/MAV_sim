import sys, os
# Add tudatpy path
if "cala" in os.getcwd():
    sys.path.append("/cala/jeremie/tudat-bundle/build/tudatpy")
else:
    sys.path.append("/mnt/c/TUDAT/tudat-bundle/build/tudatpy")
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))

from tudatpy.kernel.interface import spice
spice.load_standard_kernels()
import numpy as np

bodies = ["Sun", "Jupiter", "Venus", "Earth", "Saturn"]

for body in bodies:
    state = spice.get_body_cartesian_state_at_epoch(body, "Mars", "J2000", "None", 0)
    r = np.linalg.norm(state[0:3])
    mu = spice.get_body_gravitational_parameter(body)
    a = mu/r**2
    print(body, "%.2e m"%r, "%.2e m3/s2"%mu, "%.2e m/s2"%a, sep=", ")

phobos_r = 5989e3
phobos_mu = 10.6e15*6.6743e-11
phobos_a = phobos_mu/phobos_r**2
print("Phobos", "%.2e m"%phobos_r, "%.2e m3/s2"%phobos_mu, "%.2e m/s2"%phobos_a, sep=", ")