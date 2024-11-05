"""
Processing script for Figure 2

Reads in IVP and EVP data from the current directory and outputs fig2-processed.npy for later plotting.

Options from command line: 
    i (sys.argv[1]) index array of dispersion parameter values (alpha)
    j (sys.argv[2]) index into array of timestep sizes
"""
import numpy as np
import dedalus.public as d3
import h5py
import sys

schemes = ["sbdf1", "sbdf2", "rk222", "rk443", "sbdf3", "sbdf4"]
N = 512 
L = 10
c0 = 0.5 
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

# pick alpha and t_step to make fig of
i = int(sys.argv[1])
j = int(sys.argv[2])
alpha = alphas[i]
t_step = dts[j]

# d3 fields under complex Fourier basis to easily perform phase shifts (see uex() below)
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype = np.complex128)
xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2))
x = dist.local_grid(xbasis)
k = xbasis.wavenumbers
u_sim = dist.Field(name = 'u_sim', bases = xbasis)
u_ex = dist.Field(name = 'u_ex', bases = xbasis)
u_err = dist.Field(name = 'u_err', bases = xbasis)

# returns analytic solution 
def uex(c, alpha, t): 
    scale = c**(1/2) / (2*alpha**(1/2))
    u_ex['g'] = 3*c*np.cosh(scale*x)**(-2)
    u_ex['c'] *= np.exp(-1j*k*c*t)
    return u_ex['g']

# returns Re(lambda) from eigenvalue for a given timestep size
def eval_to_rate(eig, t_step):
    rate = np.log(np.abs(eig))/t_step
    return rate

# For loading IVP results (e.g., IVP under a real Fourier basis) 
# (NOTE: for this figure, we loaded in data that was saved at an exponentially growing cadence 
#  such that early times were sampled more frequently) 
def load_data(alphas, dts, scheme):
    data_dict = {}
    flags = np.zeros((alphas.shape[0], dts.shape[0]))
    for i, alpha in enumerate(alphas):
        data_dict[i] = {}
        for j, t_step in enumerate(dts):
            data_dict[i][j] = {}
            try:
                output_name_f = scheme + f'_output_f_{i}_{j}'
                file_in = h5py.File("figure_2_outputs_updated/" + output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
                data_dict[i][j]["ts_f"] = np.array(file_in['scales/sim_time'])
                data_dict[i][j]["u_f"] = np.array(file_in['tasks']['u_f'])
                data_dict[i][j]["f_flag"] = 1
            except:
                data_dict[i][j]["f_flag"] = 0
            flags[i, j] = data_dict[i][j]["f_flag"]
    return data_dict, flags

ts = {}
errs_ivp = {}
rates_evp = {}
flags = {} # we use flags to keep track of which (alpha, t_step) were/were not simulated

for scheme in schemes:
    # load IVP data
    data_dict, flags_out = load_data(alphas, dts, scheme)
    # load EVP data 
    print("Reading in eigenvalues assuming only the largest eigenvalues were retained")
    evals = np.load('evals_' + scheme + '_sparse_survey.npy')

    # ivp
    if not flags_out[i][j]:
        print("The selection of i", i, "j", j, "is not compatible with", scheme)
        raise
    else:
        print("Working on errors for scheme", scheme, i, j)
        ts_sim = data_dict[i][j]["ts_f"]
        ts[scheme] = ts_sim
        errs_ivp[scheme] = np.zeros(ts_sim.shape[0])
        us = data_dict[i][j]["u_f"]
        for m, t in enumerate(ts_sim):
            u_ex['g'] = uex(c0, alpha, t)
            u_sim['g'] = us[m,:]
            u_err['g'] = u_sim['g'] - u_ex['g']
            errs_ivp[scheme][m] = np.sqrt(d3.Integrate(np.abs(u_err)**2, xcoord).evaluate()['g'][0])
    # evp
    eigs = evals[i, j]
    try:
        for eig in eigs[::-1]:
            if np.logical_not(np.isclose(eig, 0+0j)):
                break
    except:
        eig = eigs
    rates_evp[scheme] = eval_to_rate(eig, t_step)


fig2_to_plot = {}
for scheme in schemes:
    fig2_to_plot[scheme] = {}

    # arranges data as iteration = 1, then 10, 20, 30, ...
    # (future users will likely want to modify this based on 
    # their output cadence strategy, which was different for each
    # scheme in our case) 
    ts_ivp = ts[scheme]
    ys_ivp = errs_ivp[scheme]
    if scheme == 'sbdf1' or scheme == 'sbdf3' or scheme == 'sbdf4':
        t_1 = ts_ivp[1]
        y_1 = ys_ivp[1]
        ts_ivp = ts_ivp[10::10]
        ys_ivp = ys_ivp[10::10]
        ts_ivp = np.concatenate(([t_1], ts_ivp))
        ys_ivp = np.concatenate(([y_1], ys_ivp))
    else:
        ts_ivp = ts_ivp[1:]
        ys_ivp = ys_ivp[1:]

    ts_evp = ts[scheme]
    lamb = rates_evp[scheme]
    # bases amplitude off of ys_ivp after one step of scheme (for us, 
    # some schemes output on iter 0 and others did not, determining
    # our truncatation above)
    ys_evp = ys_ivp[0]*np.exp(lamb*ts_evp)

    fig2_to_plot[scheme]['ivp_xaxis'] = ts_ivp/ts_ivp[-1] 
    fig2_to_plot[scheme]['ivp_yaxis'] = ys_ivp
    fig2_to_plot[scheme]['evp_xaxis'] = ts_evp/ts_evp[-1]
    fig2_to_plot[scheme]['evp_yaxis'] = ys_evp
    fig2_to_plot[scheme]['rate'] = rates_evp[scheme]

# output processed data
fig2_to_plot['(i, j)'] = (i, j)
np.save('processed_data/fig2-processed.npy', fig2_to_plot)
