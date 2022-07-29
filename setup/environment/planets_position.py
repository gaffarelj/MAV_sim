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
from datetime import datetime
from tudatpy.kernel.astro import time_conversion

bodies = ["Sun", "Jupiter", "Venus", "Earth", "Saturn"]

epoch = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime(2031, 1, 1)))
AU = 149597870700

for body in bodies:
    state = spice.get_body_cartesian_state_at_epoch(body, "Mars", "J2000", "None", epoch)
    r = np.linalg.norm(state[0:3])
    mu = spice.get_body_gravitational_parameter(body)
    a = mu/r**2
    # print(body, "%.2e m"%r, "%.2e AU"%(r/AU), "%.2e m3/s2"%mu, "%.2e m/s2"%a, sep=", ")
    print(body, "%.2e"%r, "%.2e"%(r/AU), "%.2e"%mu, "%.2e"%a, sep=" & ")

phobos_r = 5989e3
phobos_mu = 10.6e15*6.6743e-11
phobos_a = phobos_mu/phobos_r**2
# print("Phobos", "%.2e m"%phobos_r, "%.2e m3/s2"%phobos_mu, "%.2e m/s2"%phobos_a, sep=", ")
print("Phobos", "%.2e"%phobos_r, "%.2e"%(phobos_r/AU), "%.2e"%phobos_mu, "%.2e"%phobos_a, sep=" & ")