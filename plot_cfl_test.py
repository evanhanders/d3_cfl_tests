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
twidth = float(filename.split('_t')[-1].split('_')[0])
rwidth = float(filename.split('_r')[-1].split('_')[0])

data = pickle.load(open(filename,'rb'))
r_spots = data['r_spots']
θ_spots = data['θ_spots']
epsilons = data['epsilons'].transpose()

θ_plot, r_plot = np.meshgrid(r_spots, θ_spots)
#θ_plot, r_plot = np.meshgrid(θ_spots, r_spots)

print(len(r_spots), len(θ_spots))

import matplotlib.pyplot as plt
good = θ_plot <= np.pi/2
plt.pcolormesh(θ_plot, r_plot, np.log10(epsilons), cmap='viridis')
bar = plt.colorbar()
bar.set_label(r'$\log_{10}(\epsilon)$')
plt.ylabel('theta')
plt.xlabel(r'$r$')
plt.xlim(0, radius)
plt.title(r'(Lmax,Nmax) = ({}, {}) / $\sigma_{{\theta}}$, $\sigma_{{r}}$ = ({}, {}) / R = {}'.format(Lmax, Nmax, twidth, rwidth, radius))
plt.savefig('plots/growth_rates_{}.png'.format(filename.split('/')[-1].split('.pkl')[0]), dpi=300)
plt.show()
