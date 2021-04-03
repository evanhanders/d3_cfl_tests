"""
Plots data from pickle file from a cfl test

Usage:
    plot_cfl_test.py <files>...
"""
import gc
import pickle

from docopt import docopt
args   = docopt(__doc__)

import matplotlib.pyplot as plt
import numpy as np
       

plt.figure()

files = args['<files>']
for filename in files:
    radius = float(filename.split('R')[-1].split('_')[0])
    Nmax = int(filename.split('Nmax')[-1].split('.pkl')[0])
    Lmax = int(filename.split('Lmax')[-1].split('_')[0])
    twidth = float(filename.split('_t')[-1].split('_')[0])
    rwidth = float(filename.split('_r')[-1].split('_')[0])

    data = pickle.load(open(filename,'rb'))
    r_spots = data['r_spots']
    θ_spots = data['θ_spots']
    epsilons = data['epsilons'].transpose()

    #line plot
    mean_epsilon_of_r = np.mean(epsilons, axis=0)
    plt.plot(r_spots, mean_epsilon_of_r)
    plt.xlabel('r')
    plt.ylabel(r'$\epsilon$')
    plt.title(r'(Lmax,Nmax) = ({}, {}) / $\sigma_{{\theta}}$, $\sigma_{{r}}$ = ({}, {}) / R = {}'.format(Lmax, Nmax, twidth, rwidth, radius))
    plt.xlim(0, radius)

plt.savefig('plots/growth_rate_comparison.png'.format(filename.split('/')[-1].split('.pkl')[0]), dpi=300)
plt.show()
