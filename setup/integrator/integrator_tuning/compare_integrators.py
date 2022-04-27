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

# Parameters
method = "variable" # "variable" or "fixed"

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
variable_coefficients = [
    propagation_setup.integrator.heun_euler,
    propagation_setup.integrator.rkf_12,
    propagation_setup.integrator.rkf_45,
    propagation_setup.integrator.rkf_56,
    propagation_setup.integrator.rkf_78,
    propagation_setup.integrator.rkdp_87,
    propagation_setup.integrator.rkf_89,
    propagation_setup.integrator.rkv_89,
    propagation_setup.integrator.rkf_108,
    propagation_setup.integrator.rkf_1210,
    propagation_setup.integrator.rkf_1412
]

coefficients = variable_coefficients if method == "variable" else fixed_coefficients

fig = plt.figure(figsize=(10, 7))

for integ in coefficients:
    integ_name = str(integ).split(".")[-1]
    print(integ_name)

    f_evals, errors_pos, errors_vel = [], [], []
    impossible_tol = []

    # Get all the data
    res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=?", (method, integ_name))
    for row in res:
        if method == "variable" and row[4] == -1:
            if len(impossible_tol) < 2:
                impossible_tol.append(row[6])
            else:
                impossible_tol[1] = row[6]
        else:
            f_evals.append(row[3])
            errors_pos.append(row[4])
            errors_vel.append(row[5])
            print("Dt:" if method == "fixed" else "Tol:", row[7] if method == "fixed" else row[6], "Fevals:", row[3], "Error pos:", row[4], "Error vel:", row[5])
    if method == "variable" and len(impossible_tol) != 0:
        print("This tolerance range was unfeasible: [%.2e, %.2e]" % (min(impossible_tol), max(impossible_tol)))

    # Sort the lists
    try:
        f_evals, errors_pos, errors_vel = [list(tup) for tup in zip(*sorted(zip(f_evals, errors_pos, errors_vel)))]
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
plt.title("Final error in position for different %s step integrators" % method)
plt.legend()
plt.tight_layout()
plt.show()