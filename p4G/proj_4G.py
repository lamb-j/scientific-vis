import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from matplotlib.ticker import NullFormatter  # useful for `logit` scale

# Latex fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# us
error = (0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
euler_time = (21, 130, 1066, 10448, 107531, 1081077, 10493566)
rk4_time = (82, 107, 182, 310, 556, 1030, 1616)


# plot with various axes scales
plt.figure(1)

# Show how II decreases with register depth 
plt.plot(error, euler_time, 'ro', label='Euler')
plt.legend(['Euler'])
plt.plot(error, rk4_time, 'bo', label='RK4')
plt.legend(shadow=True, loc='upper left')
plt.xlim(1, 1e-6)
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Runtime (microseconds)')
plt.xlabel('Error bound')

plt.grid(True)


#plt.show()
plt.savefig("proj_4G.pdf", bbox_inches='tight')
