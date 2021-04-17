"""
Measures the growth rate due to advection of system energies using euler timestepping

Usage:
    cfl_test.py [options]

Options:
    --L=<Lmax>           The value of Lmax   [default: 14]
    --N=<Nmax>           The value of Nmax   [default: 15]
    --R=<radius>         Radius value [default: 1]

    --Nr=<N>             Number of radial points to sample CFL growth at [default: 5]
    --Ntheta=<N>         Number of latitude points to sample CFL growth at [default: 5]

    --safety=<safety>    CFL safety factor [default: 0.1]

    --theta_stdev=<ts>   Stdev of gaussian in theta direction [default: 0.04]
    --r_stdev=<ts>       Stdev of gaussian in radial direction [default: 0.04]

    --mesh=<mesh>        Processor mesh

    --shell              If flagged, use spherical shell with r_inner = 0.01R, r_outer = R
"""
import gc
import numpy as np
from mpi4py import MPI
from docopt import docopt

from dedalus.core import coords, distributor, basis, field, operators, problems, solvers, timesteppers, arithmetic
from dedalus.tools import logging
from dedalus.extras.flow_tools import CFL
from dedalus.tools.parsing import split_equation
from dedalus.extras.flow_tools import GlobalArrayReducer

import logging
logger = logging.getLogger(__name__)
args   = docopt(__doc__)

from scipy.special import erf

#Timestepper
ts = timesteppers.SBDF1

#Basis/Domains setup
Nmax = int(args['--N'])
Lmax = int(args['--L'])
radius = float(args['--R'])
dtype = np.float64
mesh = args['--mesh']
if mesh is not None:
    mesh = mesh.split(',')
    mesh = [int(mesh[0]), int(mesh[1])]
else:
    if MPI.COMM_WORLD.size > 1:
        log2 = np.log2(MPI.COMM_WORLD.size)
        if log2 == int(log2):
            mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
        logger.info("running on processor mesh={}".format(mesh))
    else:
        mesh = None

#velocity field setup
dr = radius*float(args['--r_stdev'])
dθ = np.pi*float(args['--theta_stdev'])
Nr = int(args['--Nr'])
Ntheta = int(args['--Ntheta'])
omega_magnitude = 1
vel_magnitude = 1

dealias = 3/2
c       = coords.SphericalCoordinates('φ', 'θ', 'r')
d       = distributor.Distributor((c,), mesh=mesh)
if args['--shell']:
    b       = basis.SphericalShellBasis(c, (2*(Lmax+2), Lmax+1, Nmax+1), radii=(0.01*radius, radius), dtype=dtype, dealias=(dealias, dealias, dealias))
else:
    b       = basis.BallBasis(c, (2*(Lmax+2), Lmax+1, Nmax+1), radius=radius, dtype=dtype, dealias=(dealias, dealias, dealias))
b_S2    = b.S2_basis()
φ, θ, r = b.local_grids((dealias, dealias, dealias))
φg, θg, rg = b.global_grids((dealias, dealias, dealias))

# Fields
u = field.Field(dist=d, bases=(b,), tensorsig=(c,), dtype=dtype)
a = field.Field(dist=d, bases=(b,), dtype=dtype)
omega = field.Field(dist=d, bases=(b,), dtype=dtype)
for f in [u, a, omega]:
    f.require_scales(dealias)

#For vol-avging
weight_theta = b.local_colatitude_weights(dealias)
weight_r = b.local_radial_weights(dealias)
reducer = GlobalArrayReducer(d.comm_cart)
vol_test = np.sum(weight_r*weight_theta+0*a['g'])*np.pi/(Lmax+1)/dealias
vol_test = reducer.reduce_scalar(vol_test, MPI.SUM)
vol_correction = 4*np.pi/3/vol_test

# Operations
grad = lambda A: operators.Gradient(A, c)
lap  = lambda A: operators.Laplacian(A, c)
cross = lambda A, B: arithmetic.CrossProduct(A, B)
dot = lambda A, B: arithmetic.DotProduct(A, B)
ddt = lambda A: operators.TimeDerivative(A)
LiftTau = lambda A: operators.LiftTau(A, b, -1)
radComp = lambda A: operators.RadialComponent(A)
angComp = lambda A, i=0: operators.AngularComponent(A, index=i)

lapA = lap(a)

# velocity field and initial conditions
ez = field.Field(dist=d, bases=(b,), tensorsig=(c,), dtype=dtype)
ez.set_scales(dealias)
ez['g'][1] = -np.sin(θ)
ez['g'][2] =  np.cos(θ)
ez = operators.Grid(ez).evaluate()

r_vec = field.Field(dist=d, bases=(b,), tensorsig=(c,), dtype=dtype)
r_vec.set_scales(dealias)
r_vec['g'][2] = r
r_vec = operators.Grid(r_vec).evaluate()

# Problem
problem = problems.IVP([a,])
problem.add_equation((ddt(a) + lapA - lapA, -dot(u, grad(a))))

# Solver
solver = solvers.InitialValueSolver(problem, ts, matrix_coupling=[False, False, True])
solver.stop_iteration = np.inf

reducer = GlobalArrayReducer(d.comm_cart)

noise = np.random.rand(*a['g'].shape)

def gaussian(x, mean, sigma):
    return np.exp(-(x-mean)**2/(2*sigma**2))

def run_problem(r_force, θ_force, dt, report=False): 
    omega['g'] = omega_magnitude*gaussian(r, r_force, dr)*gaussian(θ, θ_force, dθ) / np.sin(θ_force) / r_force
    u['g'] = cross(r_vec, omega*ez).evaluate()['g']


    a['g'] = 1e-1*noise
    a.require_scales(0.5)
    a['c']
    a['g']
    a.require_scales(dealias)

    E0 = reducer.reduce_scalar(np.sum(a['c']**2), MPI.SUM)
    E_now = E_prev = E0
    factor_prev = factor = 2
    epsilon_prev = epsilon = 1
    convergence = 1

    u_max = reducer.reduce_scalar(np.abs(u['g']).max(), MPI.MAX)
    logger.info('Max velocity: {:.4f}'.format(u_max))

    # Main loop
    while solver.ok and convergence > 1e-6:
        solver.step(dt)
        factor_prev = factor
        epsilon_prev = epsilon
        E_prev = E_now
        E_now = reducer.reduce_scalar(np.sum(a['c']**2), MPI.SUM)
        factor = E_now / E_prev
        epsilon = factor - 1
        convergence =  np.abs(1 - np.abs(epsilon/epsilon_prev))
        if solver.iteration % 100 == 0 and report:
            logger.info("Iter: {:05d}, Eprev: {:.03e}, Enow: {:.03e}, Epsilon: {:.03e}, convergence: {:.03e}".format(solver.iteration, E_prev, E_now, epsilon, convergence))
        a['c'] /= np.sqrt(E_now)
        E_now = 1
    return u_max, epsilon, factor

#Forcing locations
if args['--shell']:
    r_min, r_max = 0.01, 1
else:
    r_min, r_max = 5*dr/radius, 1
theta_min, theta_max = 0.5, 0.9
r_spots = radius*np.linspace(r_min, r_max, Nr)
θ_spots = np.pi * np.linspace(theta_min, theta_max, Ntheta)

#timestep
safety = float(args['--safety'])
dt = safety * radius / Lmax

#Run simulations
factors = np.zeros((len(r_spots), len(θ_spots)))
epsilons = np.zeros((len(r_spots), len(θ_spots)))
u_maxs = np.zeros((len(r_spots), len(θ_spots)))
for i, r_force in enumerate(r_spots):
    for j, θ_force in enumerate(θ_spots):
        u_max, epsilon, factor = run_problem(r_force, θ_force, dt, report=True)
        factors[i,j] = factor
        epsilons[i,j] = epsilon
        u_maxs[i,j] = u_max
        logger.info("θ/pi: {:.3f}, r/R: {:.3f}, epsilon: {:.3e}, u_max: {:.3e}".format(θ_force/np.pi, r_force/radius, epsilon, u_max))
        gc.collect()
        if not solver.ok:
            raise ValueError("Solver state not OK; something went wrong in timestepping")

#Save output
import pickle
data = {'Nmax': Nmax, 'Lmax' : Lmax, 'r_spots' : r_spots, 'θ_spots': θ_spots, 'epsilons': epsilons, 'u_maxs': u_maxs, 'factors': factors}
if args['--shell']:
    pickle.dump(data, open('dt_epsilons_shell_R{}_safety{}_stdevs_t{}_r{}_Lmax{}_Nmax{}.pkl'.format(radius, args['--safety'], args['--theta_stdev'], args['--r_stdev'], Lmax, Nmax),'wb'))
else:
    pickle.dump(data, open('dt_epsilons_ball_R{}_safety{}_stdevs_t{}_r{}_Lmax{}_Nmax{}.pkl'.format(radius, args['--safety'], args['--theta_stdev'], args['--r_stdev'], Lmax, Nmax),'wb'))

##Plot
#import matplotlib.pyplot as plt
#for i, rv in enumerate(r_spots):
#    plt.scatter(θ_spots/np.pi, epsilons[i,:], label='r_force/R = {}'.format(rv/radius))
#    plt.yscale('log')
#    plt.xlabel('theta/pi')
#    plt.ylabel(r'$\epsilons$')
#plt.title('Lmax = {} / Nmax = {} / R = {}'.format(Lmax, Nmax, radius))
#plt.legend(loc='best')
#plt.show()
#
#for i, θv in enumerate(θ_spots):
#    plt.scatter(r_spots/radius, epsilons[:, i], label=r'$\theta$_force/$\pi$ = {}'.format(θv/np.pi))
#    plt.yscale('log')
#    plt.xlabel('radius/R')
#    plt.ylabel(r'$\epsilons$')
#plt.title('Lmax = {} / Nmax = {} / R = {}'.format(Lmax, Nmax, radius))
#plt.legend(loc='best')
#plt.show()
