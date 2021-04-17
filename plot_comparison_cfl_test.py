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
from fractions import Fraction

sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=64, vmax=256))

plt.figure()

files = args['<files>']
maxv = 0
minv = np.inf
for filename in files:
    radius = float(filename.split('R')[-1].split('_')[0])
    Nmax = int(filename.split('Nmax')[-1].split('.pkl')[0])
    Lmax = int(filename.split('Lmax')[-1].split('_')[0])
    twidth = float(filename.split('_t')[-1].split('_')[0])
    rwidth = float(filename.split('_r')[-1].split('_')[0])
    safety = float(filename.split('_safety')[-1].split('_')[0])

    data = pickle.load(open(filename,'rb'))
    r_spots = data['r_spots']
    θ_spots = data['θ_spots']
    epsilons = data['epsilons'].transpose()
    print('u min/max', filename, data['u_maxs'].min(), data['u_maxs'].max())


    #line plot
    mean_epsilon_of_r = np.mean(epsilons, axis=0)
    if 'ball' in filename:
        alpha = "-1/5"
        plt.plot(r_spots, Lmax**(float(Fraction(alpha)))*mean_epsilon_of_r*r_spots**2*(safety)**(-2), color=sm.to_rgba(Nmax+1))
        plt.ylabel(r'$\epsilon$$\,\,\cdot\,\,(r\,/\,$safety)$^{2}\,$L$_{\rm{max}}$' + r'$^{{{}}}$'.format(alpha))
    else:
        plt.plot(r_spots, mean_epsilon_of_r*r_spots**2*(safety)**(-2), color=sm.to_rgba(Nmax+1))
        plt.ylabel(r'$\epsilon$$\,\,\cdot\,\,(r\,/\,$safety)$^{2}$')
    plt.xlabel('r')
    plt.xlim(r_spots.min()/2, radius)

if 'ball' in filename:
    plt.plot(r_spots, r_spots**2, label=r'$r^2$', c='r')
    plt.legend()
    plt.ylim(0, 0.4)
cbar = plt.colorbar(sm)
cbar.set_label('Nmax + 1')



#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize=8, ncol=2)
plt.savefig('plots/growth_rate_comparison.png'.format(filename.split('/')[-1].split('.pkl')[0]), dpi=300)

#plt.figure()
#
#files = args['<files>']
#for filename in files:
#    radius = float(filename.split('R')[-1].split('_')[0])
#    Nmax = int(filename.split('Nmax')[-1].split('.pkl')[0])
#    Lmax = int(filename.split('Lmax')[-1].split('_')[0])
#    twidth = float(filename.split('_t')[-1].split('_')[0])
#    rwidth = float(filename.split('_r')[-1].split('_')[0])
#
#    data = pickle.load(open(filename,'rb'))
#    r_spots = data['r_spots']
#    θ_spots = data['θ_spots']
#    epsilons = data['epsilons'].transpose()
#
#    #line plot
#    mean_epsilon_of_r = np.mean(epsilons, axis=0)
#    mean_epsilon_of_r /= mean_epsilon_of_r.max()
#    plt.plot(r_spots, mean_epsilon_of_r, label = r'(Lmax,Nmax) = ({}, {}) / $\sigma_{{\theta}}$, $\sigma_{{r}}$ = ({}, {})'.format(Lmax, Nmax, twidth, rwidth))
#    plt.xlabel('r')
#    plt.ylabel(r'$\epsilon/\epsilon_{\rm{max}}$')
#    plt.xlim(r_spots.min()/2, radius)
#
#    r_vals = np.logspace(np.log10(r_spots.min()), np.log10(r_spots.max()), 100)
#plt.plot(r_vals, 3e-4*r_vals**(-2), label=r'$r^{-2}$', c='k', lw=2)
##plt.ylim(0.1, 1.5)
#
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(fontsize=10)
#plt.savefig('plots/growth_rate_normalized_comparison.png'.format(filename.split('/')[-1].split('.pkl')[0]), dpi=300)
#plt.show()
