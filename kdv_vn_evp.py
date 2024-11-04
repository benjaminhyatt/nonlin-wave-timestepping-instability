"""
This script contains the function definitions need to solve for the 
von Neumann eigenmodes corresponding to linear perturbations of a 
KdV soliton under several IMEX schemes. 

* After using this script to perform eigenvalue problems at MULTIPLE 
resolutions (N), we implement the procedure shown in .py to decide which
eigenvalues are well-resolved. 

Following the function definitions, we also provide an example of us calling 
one of the functions to perform a parameter survey of SBDF1. 
(For higher order schemes, if all that is required is to find the fastest 
growing modes, we recommend using sparse methods--the dense methods are 
available for when the full spectrum is desired.) 

We use Dedalus** to solve each eigenvalue problem (EVP) which implements
both dense and sparse methods (both of which we include here). 

** We modify the script dedalus/core/solvers.py in Dedalus v3 in order to
   apply the desired phase shift operator S in the construction of the EVP. 
   This procedure is slightly different for multi-step vs multi-stage schemes.
   Rather than reproduce the entire Dedalus v3 build here (which
   is available at https://github.com/DedalusProject/dedalus/tree/master),
   we provide the modified solvers.py in the current repository.

Note: this script can easily be modified to perform the von Neumann analysis
for the constant background as well. To do so, you need only change the definitions
of the background fields. 

"""

import numpy as np
import dedalus.public as d3 # should point to a dedalus build with modified solvers.py


### EVPs for each scheme ###

# NOTE: for sparse solves, we supply Ne (the number of eigenmodes to be found)
# and target (a guess for where the eigenvalues (amplification factors sigma)
# are in the complex plane)

def sbdf1_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 1 # SBDF scheme
    steps = 1 
    
    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)

    problem_e = d3.EVP([u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*(u_e + alpha*t_step*dx(dx(dx(u_e)))) - (u_e - t_step*(u0*dx(u_e) + u_e*dx(u0))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]  
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]
    
    return evals, evecs

def sbdf1_dense(L, N, c, alpha, t_step):
    flag = 1 # SBDF scheme
    steps = 1 
    
    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)

    problem_e = d3.EVP([u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*(u_e + alpha*t_step*dx(dx(dx(u_e)))) - (u_e - t_step*(u0*dx(u_e) + u_e*dx(u0))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def sbdf2_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 1 # SBDF scheme
    steps = 2

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)

    problem_e = d3.EVP([u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((3/2)*u_e1 + alpha*t_step*dx(dx(dx(u_e1)))) + t_step*(2*u01*dx(u_e1)+2*dx(u01)*u_e1) - 2*u_e1 + t_step*(-u0*dx(u_e)-dx(u0)*u_e) + (1/2)*u_e = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def sbdf2_dense(L, N, c, alpha, t_step):
    flag = 1 # SBDF scheme
    steps = 2

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)

    problem_e = d3.EVP([u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((3/2)*u_e1 + alpha*t_step*dx(dx(dx(u_e1)))) + t_step*(2*u01*dx(u_e1)+2*dx(u01)*u_e1) - 2*u_e1 + t_step*(-u0*dx(u_e)-dx(u0)*u_e) + (1/2)*u_e = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def sbdf3_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 1 # SBDF scheme
    steps = 3

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-2*c*t_step))**(-2)

    problem_e = d3.EVP([u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((11/6)*u_e2 + alpha*t_step*dx(dx(dx(u_e2)))) + t_step*(3*u02*dx(u_e2)+3*dx(u02)*u_e2) - 3*u_e2 + t_step*(-3*u01*dx(u_e1)-3*dx(u01)*u_e1) + (3/2)*u_e1 + t_step*(u0*dx(u_e)+dx(u0)*u_e) - (1/3)*u_e = 0")
    problem_e.add_equation("s*u_e1 - u_e2 = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs


def sbdf3_dense(L, N, c, alpha, t_step):
    flag = 1 # SBDF scheme
    steps = 3

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-2*c*t_step))**(-2)

    problem_e = d3.EVP([u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((11/6)*u_e2 + alpha*t_step*dx(dx(dx(u_e2)))) + t_step*(3*u02*dx(u_e2)+3*dx(u02)*u_e2) - 3*u_e2 + t_step*(-3*u01*dx(u_e1)-3*dx(u01)*u_e1) + (3/2)*u_e1 + t_step*(u0*dx(u_e)+dx(u0)*u_e) - (1/3)*u_e = 0")
    problem_e.add_equation("s*u_e1 - u_e2 = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def sbdf4_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 1 # SBDF scheme
    steps = 4

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    u_e3 = dist.Field(name = 'u_e3', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-2*c*t_step))**(-2)
    u03 = dist.Field(bases = xbasis, name = 'u02')
    u03['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-3*c*t_step))**(-2)

    problem_e = d3.EVP([u_e3, u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((25/12)*u_e3 + alpha*t_step*dx(dx(dx(u_e3)))) - ((4*u_e3 - 4*t_step*(u03*dx(u_e3)+u_e3*dx(u03))) + (-3*u_e2 + 6*t_step*(u02*dx(u_e2)+u_e2*dx(u02))) + ((4/3)*u_e1 - 4*t_step*(u01*dx(u_e1)+u_e1*dx(u01))) + ((-1/4)*u_e + t_step*(u0*dx(u_e)+u_e*dx(u0)))) = 0")
    problem_e.add_equation("s*u_e2 - u_e3 = 0")
    problem_e.add_equation("s*u_e1 - u_e2 = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def sbdf4_dense(L, N, c, alpha, t_step):
    flag = 1 # SBDF scheme
    steps = 4

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    u_e3 = dist.Field(name = 'u_e3', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-2*c*t_step))**(-2)
    u03 = dist.Field(bases = xbasis, name = 'u02')
    u03['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-3*c*t_step))**(-2)

    problem_e = d3.EVP([u_e3, u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("s*((25/12)*u_e3 + alpha*t_step*dx(dx(dx(u_e3)))) - ((4*u_e3 - 4*t_step*(u03*dx(u_e3)+u_e3*dx(u03))) + (-3*u_e2 + 6*t_step*(u02*dx(u_e2)+u_e2*dx(u02))) + ((4/3)*u_e1 - 4*t_step*(u01*dx(u_e1)+u_e1*dx(u01))) + ((-1/4)*u_e + t_step*(u0*dx(u_e)+u_e*dx(u0)))) = 0")
    problem_e.add_equation("s*u_e2 - u_e3 = 0")
    problem_e.add_equation("s*u_e1 - u_e2 = 0")
    problem_e.add_equation("s*u_e - u_e1 = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, steps)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]

    return evals, evecs

def rk222_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 2 # RK scheme
    stages = 2
    gamma = (2-np.sqrt(2))/2
    delta = 1-1/(2*gamma)
    abscissae = np.array((0, gamma, 1))

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[1]*t_step))**(-2)

    problem_e = d3.EVP([u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("(u_e1 + gamma*t_step*alpha*dx(dx(dx(u_e1)))) - (u_e - gamma*t_step*(u0*dx(u_e) + u_e*dx(u0))) = 0")
    problem_e.add_equation("s*(u_e + gamma*t_step*alpha*dx(dx(dx(u_e)))) - (u_e - delta*t_step*(u0*dx(u_e) + u_e*dx(u0)) - (1-gamma)*t_step*alpha*dx(dx(dx(u_e1))) - (1-delta)*t_step*(u01*dx(u_e1) + u_e1*dx(u01))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, stages)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]
    evecs = evecs[:,np.logical_not(np.isinf(evals))]
    evals = evals[np.logical_not(np.isinf(evals))]

    return evals, evecs[N:,:]

def rk222_dense(L, N, c, alpha, t_step):
    flag = 2 # RK scheme
    stages = 2
    gamma = (2-np.sqrt(2))/2
    delta = 1-1/(2*gamma)
    abscissae = np.array((0, gamma, 1))

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[1]*t_step))**(-2)

    problem_e = d3.EVP([u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("(u_e1 + gamma*t_step*alpha*dx(dx(dx(u_e1)))) - (u_e - gamma*t_step*(u0*dx(u_e) + u_e*dx(u0))) = 0")
    problem_e.add_equation("s*(u_e + gamma*t_step*alpha*dx(dx(dx(u_e)))) - (u_e - delta*t_step*(u0*dx(u_e) + u_e*dx(u0)) - (1-gamma)*t_step*alpha*dx(dx(dx(u_e1))) - (1-delta)*t_step*(u01*dx(u_e1) + u_e1*dx(u01))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, stages)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]
    # throw out inf evals/evecs
    evecs = evecs[:,np.logical_not(np.isinf(evals))]
    evals = evals[np.logical_not(np.isinf(evals))]

    return evals, evecs[N:,:]

def rk443_sparse(L, N, c, alpha, t_step, Ne, target):
    flag = 2 # RK scheme
    stages = 4
    abscissae = np.array((0, 1/2, 2/3, 1/2, 1))

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    u_e3 = dist.Field(name = 'u_e3', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[1]*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[2]*t_step))**(-2)
    u03 = dist.Field(bases = xbasis, name = 'u03')
    u03['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[3]*t_step))**(-2)

    problem_e = d3.EVP([u_e3, u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("(u_e1 + t_step*(1/2)*alpha*dx(dx(dx(u_e1)))) - (u_e + t_step*(1/2)*(-u0*dx(u_e) - u_e*dx(u0))) = 0")
    problem_e.add_equation("(u_e2 + t_step*(1/2)*alpha*dx(dx(dx(u_e2)))) - (u_e + t_step*((11/18)*(-u0*dx(u_e) - u_e*dx(u0)) + (1/18)*(-u01*dx(u_e1) - u_e1*dx(u01))) - t_step*((1/6)*alpha*dx(dx(dx(u_e1))))) = 0")
    problem_e.add_equation("(u_e3 + t_step*(1/2)*alpha*dx(dx(dx(u_e3)))) - (u_e + t_step*((5/6)*(-u0*dx(u_e) - u_e*dx(u0)) - (5/6)*(-u01*dx(u_e1) - u_e1*dx(u01)) + (1/2)*(-u02*dx(u_e2) - u_e2*dx(u02))) - t_step*((-1/2)*alpha*dx(dx(dx(u_e1))) + (1/2)*alpha*dx(dx(dx(u_e2))))) = 0")
    problem_e.add_equation("s*(u_e + t_step*(1/2)*alpha*dx(dx(dx(u_e)))) - (u_e + t_step*((1/4)*(-u0*dx(u_e) - u_e*dx(u0)) + (7/4)*(-u01*dx(u_e1) - u_e1*dx(u01)) + (3/4)*(-u02*dx(u_e2) - u_e2*dx(u02)) - (7/4)*(-u03*dx(u_e3) - u_e3*dx(u03))) - t_step*((3/2)*alpha*dx(dx(dx(u_e1))) - (3/2)*alpha*dx(dx(dx(u_e2))) + (1/2)*alpha*dx(dx(dx(u_e3))))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_sparse(solver_e.subproblems[0], Ne, target, pars = np.array((flag, L, c, t_step, stages)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]
    evecs = evecs[:,np.logical_not(np.isinf(evals))]
    evals = evals[np.logical_not(np.isinf(evals))]

    return evals, evecs[int(3*N):,:]

def rk443_dense(L, N, c, alpha, t_step):
    flag = 2 # RK scheme
    stages = 4
    abscissae = np.array((0, 1/2, 2/3, 1/2, 1))

    dtype = np.complex128
    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype = dtype)
    xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2), dealias = 3/2)
    x = dist.local_grid(xbasis)
    dx = lambda A: d3.Differentiate(A, xcoord)

    u_e = dist.Field(name = 'u_e', bases = xbasis)
    u_e1 = dist.Field(name = 'u_e1', bases = xbasis)
    u_e2 = dist.Field(name = 'u_e2', bases = xbasis)
    u_e3 = dist.Field(name = 'u_e3', bases = xbasis)
    s = dist.Field(name = 's')
    u0 = dist.Field(bases=xbasis, name = 'u0')
    u0['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*x)**(-2)
    u01 = dist.Field(bases = xbasis, name = 'u01')
    u01['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[1]*t_step))**(-2)
    u02 = dist.Field(bases = xbasis, name = 'u02')
    u02['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[2]*t_step))**(-2)
    u03 = dist.Field(bases = xbasis, name = 'u03')
    u03['g'] = (3*c) * np.cosh((c**(1/2)/(2*alpha**(1/2)))*(x-c*abscissae[3]*t_step))**(-2)

    problem_e = d3.EVP([u_e3, u_e2, u_e1, u_e], eigenvalue = s, namespace = locals())
    problem_e.add_equation("(u_e1 + t_step*(1/2)*alpha*dx(dx(dx(u_e1)))) - (u_e + t_step*(1/2)*(-u0*dx(u_e) - u_e*dx(u0))) = 0")
    problem_e.add_equation("(u_e2 + t_step*(1/2)*alpha*dx(dx(dx(u_e2)))) - (u_e + t_step*((11/18)*(-u0*dx(u_e) - u_e*dx(u0)) + (1/18)*(-u01*dx(u_e1) - u_e1*dx(u01))) - t_step*((1/6)*alpha*dx(dx(dx(u_e1))))) = 0")
    problem_e.add_equation("(u_e3 + t_step*(1/2)*alpha*dx(dx(dx(u_e3)))) - (u_e + t_step*((5/6)*(-u0*dx(u_e) - u_e*dx(u0)) - (5/6)*(-u01*dx(u_e1) - u_e1*dx(u01)) + (1/2)*(-u02*dx(u_e2) - u_e2*dx(u02))) - t_step*((-1/2)*alpha*dx(dx(dx(u_e1))) + (1/2)*alpha*dx(dx(dx(u_e2))))) = 0")
    problem_e.add_equation("s*(u_e + t_step*(1/2)*alpha*dx(dx(dx(u_e)))) - (u_e + t_step*((1/4)*(-u0*dx(u_e) - u_e*dx(u0)) + (7/4)*(-u01*dx(u_e1) - u_e1*dx(u01)) + (3/4)*(-u02*dx(u_e2) - u_e2*dx(u02)) - (7/4)*(-u03*dx(u_e3) - u_e3*dx(u03))) - t_step*((3/2)*alpha*dx(dx(dx(u_e1))) - (3/2)*alpha*dx(dx(dx(u_e2))) + (1/2)*alpha*dx(dx(dx(u_e3))))) = 0")
    solver_e = problem_e.build_solver(ncc_cutoff = 1e-16)
    solver_e.solve_dense(solver_e.subproblems[0], pars = np.array((flag, L, c, t_step, stages)))
    idx = np.abs(solver_e.eigenvalues).argsort()[::-1]
    evals = solver_e.eigenvalues[idx]
    evecs = solver_e.eigenvectors[:,idx]
    evecs = evecs[:,np.logical_not(np.isinf(evals))]
    evals = evals[np.logical_not(np.isinf(evals))]

    return evals, evecs[int(3*N):,:]

### Example usage ###

# Fixed parameters
c = 0.5
L = 10

# Parameters to survey
alphas = np.logspace(np.log10(3e-3), np.log10(2e-2), 10)
dts = np.logspace(np.log10(6e-4), np.log10(5e-2), 22)

# Resolution
N = 512

dtype = np.complex128
evals_sbdf1 = np.zeros((alphas.shape[0], dts.shape[0], N), dtype=dtype)
evecs_sbdf1 = np.zeros((alphas.shape[0], dts.shape[0], N, N), dtype=dtype)

for i, alpha in enumerate(alphas):
    for j, t_step in enumerate(dts):
        print("Now solving with", i, f"alpha = {alphas[i]}", j, f"t_step = {dts[j]}")
        evals_sbdf1[i, j, :], evecs_sbdf1[i, j, :, :] = sbdf1_dense(L, N, c, alpha, t_step)

# Outputs
print("Outputting all eigenvalues and all eigenvectors--")
print("If you only want to save those corresponding to the fastest growing modes, uncomment next two lines")
# evals_sbdf1 = evals_sbdf1[:, :, 0]
# evecs_sbdf1 = evecs_sbdf1[:, :, :, 0]
np.save("evals_sbdf1_dense_survey.npy", evals_sbdf1)
np.save("evecs_sbdf1_dense_survey.npy", evecs_sbdf1)
