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

# Standard imports
import numpy as np
import sqlite3
from matplotlib import pyplot as plt

# Tudatpy imports
from tudatpy.kernel.numerical_simulation import propagation_setup

# Connect to the database
con = sqlite3.connect(sys.path[0]+"/setup/integrator/integrator_tuning/database.db")
cur = con.cursor()

# List of integrators
fixed_coefficients = [
    propagation_setup.integrator.euler_forward,
    propagation_setup.integrator.rk_4,
    propagation_setup.integrator.explicit_mid_point,
    propagation_setup.integrator.explicit_trapezoid_rule,
    propagation_setup.integrator.ralston,
    propagation_setup.integrator.rk_3,
    propagation_setup.integrator.ralston_3,
    propagation_setup.integrator.SSPRK3,
    propagation_setup.integrator.ralston_4,
    propagation_setup.integrator.three_eight_rule_rk_4
]

fig = plt.figure(figsize=(10, 7))

for integ in fixed_coefficients:
    integ_name = str(integ).split(".")[-1]
    print(integ_name)

    dts, f_evals, errors_pos, errors_vel = [], [], [], []

    # Get all the data
    res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=?", ("fixed",integ_name))
    for row in res:
        dts.append(row[7])
        f_evals.append(row[3])
        errors_pos.append(row[4])
        errors_vel.append(row[5])
        print("Dt:", row[7], "Fevals:", row[3], "Error pos:", row[4], "Error vel:", row[5])

    # Sort the lists
    try:
        dts, f_evals, errors_pos, errors_vel = [list(tup) for tup in zip(*sorted(zip(dts, f_evals, errors_pos, errors_vel)))]
    except ValueError:
        print("No data for this integrator")
        continue

    # Plot the data
    plt.plot(f_evals, errors_pos, label=integ_name, marker="o", linestyle="-")


# Close database connection and show plot
con.close()
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of function evaluations [-]")
plt.ylabel("Final error in position [m]")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()