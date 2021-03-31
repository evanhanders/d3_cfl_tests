"""
Plots data from pickle file from a cfl test

Usage:
    plot_cfl_test.py <file>
"""
from docopt import docopt
args   = docopt(__doc__)


import gc
import numpy as np
       
import pickle

filename = args['<file>']
radius = float(filename.split('R')[-1].split('_')[0])
Nmax = int(filename.split('Nmax')[-1].split('.pkl')[0])
Lmax = int(filename.split('Lmax')[-1].split('_')[0])

data = pickle.load(open(filename,'rb'))
r_spots = data['r_spots']
θ_spots = data['θ_spots']
epsilons = data['epsilons']

print(len(r_spots), len(θ_spots))

import matplotlib.pyplot as plt
plt.scatter(θ_spots/np.pi, epsilons[0,:], label='r_force/R = {}'.format(r_spots[0]/radius))
plt.scatter(θ_spots/np.pi, epsilons[-1,:], label='r_force/R = {}'.format(r_spots[-1]/radius))
plt.yscale('log')
plt.xlabel('theta/pi')
plt.ylabel(r'$\epsilon$')
plt.title('Lmax = {} / Nmax = {} / R = {}'.format(Lmax, Nmax, radius))
plt.legend(loc='best')
plt.show()

for i, θv in enumerate(θ_spots):
    plt.scatter(r_spots/radius, epsilons[:, i], label=r'$\theta$_force/$\pi$ = {}'.format(θv/np.pi))
    plt.yscale('log')
    plt.xlabel('radius/R')
    plt.ylabel(r'$\epsilon$')
plt.title('Lmax = {} / Nmax = {} / R = {}'.format(Lmax, Nmax, radius))
#plt.legend(loc='best')
plt.show()
