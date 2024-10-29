"""
This is a script that uses Dedalus to simulate soliton solutions of
the nonlinear KdV equation. 

Options from command line:
    i (sys.argv[1]) index array of dispersion parameter values (alpha)
    j (sys.argv[2]) index into array of timestep sizes

Options in script:
    L (line 54)           size of 1D domain
    c0 (line 55)          soliton parameter c for initial condition
    N (line 56)           number of basis polynomials
    timestepper (line 59) choice of IMEX timestepping scheme     
    T_end (line 79)       specified stop time of simulation 
    basis (line 84)       choice of spectral basis
    out_folder (line 86)  name of folder to output analysis tasks

Outputs:
    u            snapshots of the scalar velocity field in grid space
    int_u        snapshots of the integral of the scalar velocity field
    int_u_sq     snapshots of the L2 norm of the scalar velocity field   

To specify other outputs or alter the output cadence strategy, we refer
the reader to the Dedalus documentation which can be found at:
https://dedalus-project.readthedocs.io/en/latest/
"""

# import dependencies 
import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)
import sys

### Specify options ###

# parameter ranges tested in this project
alphas = np.logspace(np.log10(3e-3), np.log10(2e-2), 10) 
dts = np.logspace(np.log10(6e-4), np.log10(5e-2), 22) 
# IMEX schemes considered in this work
schemes = ['sbdf1', 'sbdf2', 'sbdf3', 'sbdf4', 'rk222', 'rk443']
# spectral bases considered in this work
bases = ['fourier', 'chebyshev']

# select parameters to simulate
i = int(sys.argv[1])
j = int(sys.argv[2])

print("Read i =", i, " and j =", j)

alpha = alphas[i]
t_step = dts[j]

L = 10
c0 = 0.5
N = 512

s = 1 # (e.g., for sbdf2)
scheme = schemes[s]

# Here we specify expected blowup/decay times for some of 
# the schemes that we looked at with an asymptotic analysis
# in this work--for schemes that have not been so analyzed, 
#one should specify their own desired stop time
def T_late(scheme, c0, alpha, t_step):
    if scheme == 'sbdf1':
        T = (35/34) * c0**(-3) * alpha * t_step**(-1)
    elif scheme == 'sbdf2':
        T = (35/86) * c0**(-6) * alpha**2 * t_step**(-3)
    elif scheme == 'rk222':
        T = (5005/(355045-245436*np.sqrt(2))) * c0**(-6) * alpha**2 * t_step**(-3)
    elif scheme == 'rk443':
        T = 3 * (30030/77069) * c0**(-6) * alpha**2 * t_step**(-3)
    else:
        print("No T_end has been specified for scheme", scheme)
        raise
    return T

T_end = 3 * T_late(scheme, c0, alpha, t_step)

print("alpha = " + f"{alpha:.5f}, " + "t_step = " + f"{t_step:.6f}, " + "Enforced stop time = " + f"{T_end:.2f}")

b = 0 # (e.g., for real fourier)
basis = bases[b]
basis_char = basis[0] 
out_folder = scheme + '_output_' + basis_char + f'_{i}_{j}'

### Begin Dedalus setup ###
dealias = 3/2 
dtype = np.float64
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=dtype)

if basis == 'fourier':
    # specify real fourier basis
    xbasis = d3.RealFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = dealias)
    dx = lambda A: d3.Differentiate(A, xcoord)
    x = dist.local_grid(xbasis)
    u = dist.Field(name = 'u', bases = xbasis)
    ux = dx(u)
    uxx = dx(ux)
    # set up IVP
    problem = d3.IVP([u], namespace = locals())
    problem.add_equation("dt(u) + alpha*dx(uxx) = -u*ux")

elif basis == 'chebyshev':
    # specify chebyshev basis
    xbasis = d3.ChebyshevT(xcoord, size = N, bounds = (-L/2, L/2), dealias = dealias)
    dx = lambda A: d3.Differentiate(A, xcoord)
    x = dist.local_grid(xbasis)
    lift_basis = xbasis.derivative_basis()
    lift = lambda A: d3.Lift(A, lift_basis, -1) 
    u = dist.Field(name = 'u', bases = xbasis)
    tau1 = dist.Field(name = 'tau1')
    tau2 = dist.Field(name = 'tau2')
    tau3 = dist.Field(name = 'tau3')
    ux = dx(u) + lift(tau1)
    uxx = dx(ux) + lift(tau2)
    # set up IVP
    problem = d3.IVP([u, tau1, tau2, tau3], namespace = locals())
    problem.add_equation("dt(u) + alpha*dx(uxx) + lift(tau3) = -u*ux")
    problem.add_equation("u(x = -L/2) - u(x = L/2) = 0")
    problem.add_equation("ux(x = -L/2) - ux(x = L/2) = 0")
    problem.add_equation("uxx(x = -L/2) - uxx(x = L/2) = 0") 
else:
    print(basis, "is not implemented")
    raise 

# specify timestepper
if scheme == 'sbdf1':
    timestepper = d3.SBDF1
elif scheme == 'sbdf2':
    timestepper = d3.SBDF2
elif scheme == 'sbdf3':
    timestepper = d3.SBDF3
elif scheme == 'sbdf4':
    timestepper = d3.SBDF4
elif scheme == 'rk222':
    timestepper = d3.RK222
elif scheme == 'rk443':
    timestepper = d3.RK443
else:
    print(scheme, "not implemented")
    raise

# build IVP solver
solver = problem.build_solver(timestepper)
solver.sim_time = 0
solver.stop_sim_time = T_end

# specify initial condition
u['g'] = (3*c0) * np.cosh( np.sqrt(c0/(4*alpha)) * x) ** (-2)

# setup analysis tasks
analysis = solver.evaluator.add_file_handler(out_folder, iter = 1000)
analysis.add_task(u, layout='g', name='u')
analysis.add_task(d3.Integrate(u, 'x'), name='int_u')
analysis.add_task(d3.Integrate(u**2, 'x'), name = 'int_u_sq')

while solver.proceed:
    if solver.iteration % 100 == 0:
        logger.info('Iteration=%i, Time=%e' %(solver.iteration, solver.sim_time))
    solver.step(t_step)
    # check for blow-up
    if not np.isfinite(u['g']).all():
        solver.evaluate_handlers()
        logger.info('Iteration=%i, Time=%e, solver encountered non-finite value on grid - breaking loop' %(solver.iteration, solver.sim_time))
        break

