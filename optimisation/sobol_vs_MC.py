import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc


x_lims, y_lims = [-1, 1], [0, 10]

# Generate a random sample of size n
n = 150
seed = 123987
np.random.seed(seed)
# get log base 2 of n rounded to higher integer
log_n = int(np.ceil(np.log2(1*n)))
x_MC = np.random.uniform(x_lims[0], x_lims[1], int(1*n))
y_MC = np.random.uniform(y_lims[0], y_lims[1], int(1*n))
print("Generating %i samples from MC" % int(1*n))

# Generate samples from Sobol sequence
sampler = qmc.Sobol(d=2, scramble=True, seed=seed)
samples = sampler.random_base2(m=log_n) # m=3 means 2^3=8 samples between 0 and 1
# Make values fit in limits
x_sobol = x_lims[0] + samples[:, 0] * (x_lims[1] - x_lims[0])
y_sobol = y_lims[0] + samples[:, 1] * (y_lims[1] - y_lims[0])
print("Generating %i samples from Sobol" % int(2**log_n))

# Plot the sample
plt.scatter(x_MC, y_MC, label='Monte Carlo')
plt.scatter(x_sobol, y_sobol, label='Sobol')
plt.xlim(x_lims)
plt.ylim(y_lims)
plt.grid()
plt.legend()
plt.xlabel('x'), plt.ylabel('y')
plt.tight_layout()
plt.show()
