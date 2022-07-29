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
method = "fixed" # "variable" or "fixed"
allowed_errors = [5e3, 5]

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

coefficients = variable_coefficients if method == "variable" else fixed_coefficients#+variable_coefficients

fig1 = plt.figure(figsize=(10, 6))
fig2 = plt.figure(figsize=(10, 6))

i_plotted = 0
for integ in coefficients:
    integ_name = str(integ).split(".")[-1]
    print(integ_name)

    f_evals, errors_pos, errors_vel = [], [], []
    impossible_tol = []

    # Get all the data
    res = cur.execute("SELECT * FROM integrator WHERE method=? AND coefficients=? ORDER BY f_evals ASC", (method, integ_name))
    for row in res:
        if method == "variable" and row[4] == -1:
            if len(impossible_tol) < 2:
                impossible_tol.append(row[6])
            else:
                impossible_tol[1] = row[6]
        else:
            if row[4] != 999999999 and row[5] != 999999999:
                if (row[4] >= 50 and row[5] >= 0.05 and row[3] < 3e4) or not method == "variable":
                    f_evals.append(row[3])
                    errors_pos.append(row[4])
                    errors_vel.append(row[5])
                    i_plotted += 1
                    print("Dt:" if method == "fixed" else "Tol:", row[7] if method == "fixed" else row[6], "Fevals:", row[3], "Error pos:", row[4], "Error vel:", row[5])
            elif row[4] == 999999999:
                print("Impossible tol:", row[6])
    if method == "variable" and len(impossible_tol) != 0:
        print("This tolerance range was unfeasible: [%.2e, %.2e]" % (min(impossible_tol), max(impossible_tol)))

    # Sort the lists
    try:
        f_evals, errors_pos, errors_vel = [list(tup) for tup in zip(*sorted(zip(f_evals, errors_pos, errors_vel)))]
    except ValueError:
        print("No data for this integrator")
        continue

    # Plot the data
    plt.figure(fig1)
    plt.plot(f_evals, errors_pos, label=integ_name, marker="o", linestyle="-")
    plt.figure(fig2)
    plt.plot(f_evals, errors_vel, label=integ_name, marker="o", linestyle="-")

print("Total of %i points plotted" % i_plotted)

# Close database connection
con.close()

# Finish plot
plt.figure(fig1)
plt.xscale("log"), plt.yscale("log")
xlim, ylim = plt.xlim(), plt.ylim()
plt.hlines(y=allowed_errors[0], xmin=xlim[0]/1e6, xmax=xlim[1]*1e6, color="orange", linestyle="--")
plt.xlim(xlim[0], xlim[1])
plt.xlabel("Number of function evaluations [-]")
plt.ylabel("Final error in position [m]")
plt.grid()
plt.title("Final error in position for different %s step integrators" % method)
plt.legend()
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/integrator/compare_integrators_%s_pos.pdf" % method)

plt.figure(fig2)
plt.xscale("log"), plt.yscale("log")
xlim, ylim = plt.xlim(), plt.ylim()
plt.hlines(y=allowed_errors[1], xmin=xlim[0]/1e6, xmax=xlim[1]*1e6, color="orange", linestyle="--")
plt.xlim(xlim[0], xlim[1])
plt.xlabel("Number of function evaluations [-]")
plt.ylabel("Final error in velocity [m/a]")
plt.grid()
plt.title("Final error in velocity for different %s step integrators" % method)
plt.legend()
plt.tight_layout()
plt.savefig(sys.path[0]+"/plots/setup/integrator/compare_integrators_%s_vel.pdf" % method)
plt.close()