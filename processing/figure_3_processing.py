"""
Processing script for Figure 3

Reads in IVP and EVP data from the current directory and outputs fig3-processed.npy for later plotting.
"""
import numpy as np
import dedalus.public as d3
import h5py

schemes = ["sbdf3", "sbdf4"]
N = 512
L = 10
c0 = 0.5
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10)
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22)

# d3 fields in ComplexFourier, in order to cleanly calculate exact solution
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

# obtain an exponential fit of the error growth at early times
def exp_fit_err(ts, ys):
    tol = 1e-3
    flag = True

    # establish a high value around/below where the error is expected to saturate
    hi = np.min((0.5, np.max(ys[np.isfinite(ys)])))
    # establish low value to start fitting from
    lo = np.min(ys[np.isfinite(ys)][2:])

    # fit over this first interval
    ts_0 = ts[np.logical_and(ys <= hi, ys >= lo)]
    ys_0 = ys[np.logical_and(ys <= hi, ys >= lo)]
    p, res, rank, sv, rcond = np.polyfit(ts_0, np.log(ys_0), 1, full = True)

    # if applicable, refine the fit
    while flag and (res > tol*ts_0.shape[0]):
        lo *= 1.5 
        if lo > 1e-1*hi:
            break

        ts_0 = ts[np.logical_and(ys <= hi, ys >= lo)]
        ys_0 = ys[np.logical_and(ys <= hi, ys >= lo)]
        p, res, rank, sv, rcond = np.polyfit(ts_0, np.log(ys_0), 1, full = True)

    return p[0]

# For loading IVP results (e.g., IVP under a real Fourier basis) 
# (NOTE: for this figure, we were able to use a fixed output cadence since
# SBDF3 and SBDF4 both terminate in a relatively short time frame)
def load_data(alphas, dts, scheme):
    data_dict = {}
    flags = np.zeros((alphas.shape[0], dts.shape[0]))
    for i, alpha in enumerate(alphas):
        data_dict[i] = {}
        for j, t_step in enumerate(dts):
            data_dict[i][j] = {}
            try:
                output_name_f = scheme + f'_output_f_{i}_{j}'
                file_in = h5py.File(scheme + "_outputs_N_512_f/" + output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
                data_dict[i][j]["ts_f"] = np.array(file_in['scales/sim_time'])
                data_dict[i][j]["u_f"] = np.array(file_in['tasks']['u_f'])
                data_dict[i][j]["f_flag"] = 1
            except:
                data_dict[i][j]["f_flag"] = 0
            flags[i, j] = data_dict[i][j]["f_flag"]
    return data_dict, flags

errs = {}
rates_sim = {}
rates_evp = {}
flags = {}

for scheme in schemes:

    errs[scheme] = {}
    rates_sim[scheme] = np.zeros((alphas.shape[0], dts.shape[0]))
    rates_evp[scheme] = np.zeros((alphas.shape[0], dts.shape[0]))
    flags[scheme] = np.zeros((alphas.shape[0], dts.shape[0]))
    
    # Load IVP data
    data_dict, flags_out = load_data(alphas, dts, scheme)
    flags[scheme] = flags_out
    
    # Load EVP data
    print("Reading in eigenvalues assuming only the largest eigenvalues were retained")
    evals = np.load('evals_' + scheme + '_dense_survey.npy')

    # Calculate errors and fit rates, and also convert eigenvalues to rates
    for i, alpha in enumerate(alphas):
        errs[scheme][i] = {}
        for j, t_step in enumerate(dts):
            # IVP
            if data_dict[i][j]["f_flag"]:
                print(i, j)
                ts_sim = data_dict[i][j]["ts_f"]
                errs[scheme][i][j] = np.zeros(ts_sim.shape[0])
                us = data_dict[i][j]["u_f"]
                for m, t in enumerate(ts_sim):
                    u_ex['g'] = uex(c0, alpha, t)
                    u_sim['g'] = us[m,:]
                    u_err['g'] = u_sim['g'] - u_ex['g']
                    errs[scheme][i][j][m] = np.sqrt(d3.Integrate(u_err**2, xcoord).evaluate()['g'][0])
                # fit exp profile with refinement
                rates_sim[scheme][i, j] = exp_fit_err(ts_sim, errs[scheme][i][j])
            # EVP
            eigs = evals[i, j]
            try:
                for eig in eigs[::-1]:
                    if np.logical_not(np.isclose(eig, 0+0j)):
                        break
            except:
                eig = eigs
            rates_evp[scheme][i, j] = eval_to_rate(eig, t_step)

    # Removing if not well-resolved
    flags[scheme][0, -1] = 0
    flags[scheme][1, -1] = 0
    
fig3_to_plot = {}
for scheme in schemes:
    fig3_to_plot[scheme] = {}

    min_y = 1

    # IVP
    for i, alpha in enumerate(alphas):
        dts_plot = []
        rates_ivp_plot = []
        for j, t_step in enumerate(dts):
            if flags[scheme][i][j]:
                rates_ivp_plot.append(rates_sim[scheme][i][j])
                if rates_sim[scheme][i][j] < min_y:
                    min_y = rates_sim[scheme][i][j]
                dts_plot.append(t_step)

        fig3_to_plot[scheme]['IVP' + str(i)] = {}
        fig3_to_plot[scheme]['IVP' + str(i)]['IVP_xaxis'] = dts_plot
        fig3_to_plot[scheme]['IVP' + str(i)]['IVP_yaxis'] = rates_ivp_plot
    
    # EVP
    for i, alpha in reversed(list(enumerate(alphas))):
        dts_plot = []
        rates_evp_plot = []
        for j, t_step in enumerate(dts):
            if flags[scheme][i][j]:
                rates_evp_plot.append(rates_evp[scheme][i][j])
                if rates_evp[scheme][i][j] < min_y:
                    min_y = rates_evp[scheme][i][j]
                dts_plot.append(t_step)

        fig3_to_plot[scheme]['EVP' + str(i)] = {}
        fig3_to_plot[scheme]['EVP' + str(i)]['EVP_xaxis'] = dts_plot
        fig3_to_plot[scheme]['EVP' + str(i)]['EVP_yaxis'] = rates_evp_plot

    # Trend line
    y = np.array(dts_plot)
    yinit = rates_evp_plot[0] / 10
    idxlow = 0
    idxhi = np.where(((yinit/y[0]**(-1)) * y**(-1)) < min_y)[0][1]
    x_midpt_tl = np.logspace(np.log10((dts_plot[idxlow:idxhi])[0]), np.log10((dts_plot[idxlow:idxhi])[-1]), 3)[1]
    y_midpt_tl = np.logspace(np.log10((((yinit/y[0]**(-1)) * y**(-1))[idxlow:idxhi])[0]), np.log10((((yinit/y[0]**(-1)) * y**(-1))[idxlow:idxhi])[-1]), 3)[1]
    
    fig3_to_plot[scheme]['trend_y'] = y
    fig3_to_plot[scheme]['trend_misc'] = (yinit, idxlow, idxhi, x_midpt_tl, y_midpt_tl)


# output processed data
np.save('processed_data/fig3-processed.npy', fig3_to_plot)
