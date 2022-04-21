
"""
Design variables (17-21):
 - (1) Launch angle of SRM 1
 - (1) Launch angle of SRM 2
 - (3-7) SRM geometry parameters of stage 1 (to be defined before optimisation)
 - (2) SRM geometry parameters of stage 2 (spherical)
 - (10) TVC deflection for 5 nodes over burn time in both planes

Design variables of each SRM geometric model:

* Tubular:
    - L: length [m] ~ [0.3, 1.25]
    - R_o: outer radius [m] ~ [0.1, 0.285]
    - R_i_frac: inner radius = R_i_frac * R_o [-] ~ [0.2, 0.9]

* Rod and tube:
    - L: length [m] ~ [0.3, 1.25]
    - R_o: outer radius [m] ~ [0.1, 0.285]
    - R_mid_frac: intermediate radius = R_mid_frac * R_o [-] ~ [0.2, 0.9]
    - R_i_frac: inner radius = R_i_frac * R_mid [-] ~ [0.2, 0.9]

* Spherical:
    - R_o_frac: outer radius = R_o_frac * R_o_stage1 [m] ~ [0.3, 1.0]
    - R_i_frac: inner radius = R_i_frac * R_o [-] ~ [0.15, 0.8]

* Multi-fin:
    - L: length [m] ~ [0.3, 1.25]
    - R_o: outer radius [m] ~ [0.1, 0.285]
    - R_i_frac: inner radius = R_i_frac * R_o [-] ~ [0.2, 0.9]
    - N_f: number of fins [-] ~ [3-15]
    - L_f_frac: fin length = L_f_frac * R_i [-] ~ [0.25, 0.75]
    - w_f_frac: fin thickness = w_f_frac * 2*pi*(R_i-L_f)/N_f [m] ~ [0.35, 0.9]

* Anchor:
    - TODO
"""

import sys
# Set path to uppermost project level
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "MAV_sim":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import numpy as np
from matplotlib import pyplot as plt

from thrust.models.multi_fin import multi_fin_SRM
from thrust.models.spherical import spherical_SRM
from thrust.models.tubular import tubular_SRM
from thrust.models.rod_and_tube import rod_and_tube_SRM

MC_runs = 0
fact_run = True
test_multi_fin_dep_vars = True
test_spherical_dep_vars = True
test_tubular_dep_vars = True
test_rod_and_tube_dep_vars = True

operators_to_run = {
    "min": min,
    "0.25": lambda arr: arr[0]+0.25*(arr[1]-arr[0]),
    "0.5": lambda arr: arr[0]+0.5*(arr[1]-arr[0]),
    "0.75": lambda arr: arr[0]+0.75*(arr[1]-arr[0]),
    "max": max
}

if test_multi_fin_dep_vars:
    des_vars_range = [
        [0.3, 0.1, 0.2, 3, 0.25, 0.35],
        [1.25, 0.285, 0.9, 15, 0.75, 0.9]
    ]
    base_vals = [1, 0.15, 0.5, 7, 0.5, 0.5]
    des_vars_names = ["L", "R_o", "R_{i_{frac}}", "N_f", "L_{f_{frac}}", "w_{f_{frac}}"]
    des_vars_type = [float, float, float, int, float, float]
    if fact_run:
        fig, axes = plt.subplots(len(des_vars_range[0])-1, len(operators_to_run.keys()), figsize=(15, 15), subplot_kw=dict(projection='polar'))
        for i, des_var_range in enumerate(np.asarray(des_vars_range).T):
            if i != 0: # skip SRM length (can be whatever)
                print(i, des_var_range)
                for j, (op_name, op) in enumerate(operators_to_run.items()):

                    des_var = base_vals.copy()
                    des_var[i] = des_vars_type[i](op(des_var_range))
                    print(op_name, des_var)

                    L, R_o, R_i_frac, N_f, L_f_frac, w_f_frac = des_var
                    R_i = R_i_frac * R_o
                    L_f = L_f_frac * R_i
                    w_f = 2*np.pi*(R_i-L_f)/N_f*w_f_frac

                    ax = axes[i-1, j]
                    multi_fin_SRM(R_o, R_i, N_f, w_f, L_f, L, False).plot_geometry(ax_in=ax)
                    ax.set_title("$%s$ = %.2f" % (des_vars_names[i], des_var[i]), fontsize=14)
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/optimisation/des_vars_verif/factorial/multi_fin_SRM.pdf")
        plt.close()
    
    for i in range(MC_runs):
        rnd_des_var = []
        for j in range(len(des_vars_range[0])):
            if des_vars_type[j] == int:
                rnd_des_var.append(np.random.randint(des_vars_range[0][j], des_vars_range[1][j]))
            else:
                rnd_des_var.append(np.random.uniform(des_vars_range[0][j], des_vars_range[1][j]))

        L, R_o, R_i_frac, N_f, L_f_frac, w_f_frac = base_vals
        R_i = R_i_frac * R_o
        L_f = L_f_frac * R_i
        w_f = 2*np.pi*(R_i-L_f)/N_f*w_f_frac
        geo_vars = [R_o, R_i, N_f, w_f, L_f, L]

        multi_fin_SRM(*geo_vars, False).plot_geometry()
        plt.savefig(sys.path[0]+"/optimisation/des_vars_verif/MC/multi_fin_SRM_"+"_".join(["%.2f" % val for val in rnd_des_var])+".pdf")
        plt.close()

if test_spherical_dep_vars:
    des_vars_range = [
        [0.1, 0.15],
        [1.0, 0.8]
    ]
    base_vals = [1.0, 0.35]
    des_vars_names = ["R_o", "R_{i_{frac}}"]
    des_vars_type = [float, float]
    if fact_run:
        fig, axes = plt.subplots(len(des_vars_range[0])-1, len(operators_to_run.keys()), figsize=(15, 4), subplot_kw=dict(projection='polar'))
        for i, des_var_range in enumerate(np.asarray(des_vars_range).T):
            if i != 0: # skip SRM length (can be whatever)
                print(i, des_var_range)
                for j, (op_name, op) in enumerate(operators_to_run.items()):

                    des_var = base_vals.copy()
                    des_var[i] = des_vars_type[i](op(des_var_range))
                    print(op_name, des_var)

                    R_o, R_i_frac = des_var
                    R_i = R_i_frac * R_o

                    # ax = axes[i-1, j]
                    ax = axes[j]
                    spherical_SRM(R_o, R_i).plot_geometry(ax_in=ax)
                    ax.set_title("$%s$ = %.2f" % (des_vars_names[i], des_var[i]), fontsize=14)
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/optimisation/des_vars_verif/factorial/spherical_SRM.pdf")
        plt.close()

if test_tubular_dep_vars:
    des_vars_range = [
        [0.3, 0.05, 0.2],
        [1.25, 0.285, 0.9]
    ]
    base_vals = [1.0, 0.15, 0.5]
    des_vars_names = ["L", "R_o", "R_{i_{frac}}"]
    des_vars_type = [float, float, float]
    if fact_run:
        fig, axes = plt.subplots(len(des_vars_range[0])-1, len(operators_to_run.keys()), figsize=(15, 7), subplot_kw=dict(projection='polar'))
        for i, des_var_range in enumerate(np.asarray(des_vars_range).T):
            if i != 0: # skip SRM length (can be whatever)
                print(i, des_var_range)
                for j, (op_name, op) in enumerate(operators_to_run.items()):

                    des_var = base_vals.copy()
                    des_var[i] = des_vars_type[i](op(des_var_range))
                    print(op_name, des_var)

                    L, R_o, R_i_frac = des_var
                    R_i = R_i_frac * R_o

                    ax = axes[i-1, j]
                    tubular_SRM(R_o, R_i, L).plot_geometry(ax_in=ax)
                    ax.set_title("$%s$ = %.2f" % (des_vars_names[i], des_var[i]), fontsize=14)
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/optimisation/des_vars_verif/factorial/tubular_SRM.pdf")
        plt.close()

if test_rod_and_tube_dep_vars:
    des_vars_range = [
        [0.3, 0.05, 0.2, 0.2],
        [1.25, 0.285, 0.9, 0.9]
    ]
    base_vals = [1.0, 0.15, 0.5, 0.5]
    des_vars_names = ["L", "R_o", "R_{mid_{frac}}", "R_{i_{frac}}"]
    des_vars_type = [float, float, float, float]
    if fact_run:
        fig, axes = plt.subplots(len(des_vars_range[0])-1, len(operators_to_run.keys()), figsize=(15, 10), subplot_kw=dict(projection='polar'))
        for i, des_var_range in enumerate(np.asarray(des_vars_range).T):
            if i != 0: # skip SRM length (can be whatever)
                print(i, des_var_range)
                for j, (op_name, op) in enumerate(operators_to_run.items()):

                    des_var = base_vals.copy()
                    des_var[i] = des_vars_type[i](op(des_var_range))
                    print(op_name, des_var)

                    L, R_o, R_mid_frac, R_i_frac = des_var
                    R_mid = R_mid_frac * R_o
                    R_i = R_i_frac * R_mid

                    ax = axes[i-1, j]
                    rod_and_tube_SRM(R_o, R_mid, R_i, L).plot_geometry(ax_in=ax)
                    ax.set_title("$%s$ = %.2f" % (des_vars_names[i], des_var[i]), fontsize=14)
        plt.tight_layout()
        plt.savefig(sys.path[0]+"/optimisation/des_vars_verif/factorial/rod_and_tube_SRM.pdf")
        plt.close()