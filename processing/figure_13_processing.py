""" 
Processing script for Figure 13

Reads in VN EVP data from current directory.
Outputs fig13-processed.npy for later plotting.
"""
import numpy as np
import dedalus.public as d3
import h5py

schemes = ['sbdf1', 'sbdf2', 'rk222', 'rk443']

N = 512 
L = 10
c0 = 0.5 
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

t_comp = 10*np.max(dts)

# d3 fields in ComplexFourier, in order to cleanly calculate exact solution
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype = np.complex128)
xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2))
x = dist.local_grid(xbasis)
dx = lambda A: d3.Differentiate(A, xcoord)
k = xbasis.wavenumbers
u_sim = dist.Field(name = 'u_sim', bases = xbasis)
u_ex = dist.Field(name = 'u_ex', bases = xbasis)
u_err = dist.Field(name = 'u_err', bases = xbasis)
u_err_NL = dist.Field(name = 'u_err_NL', bases = xbasis)

# returns exact solution
def uex(c, alpha, t): 
    scale = c**(1/2) / (2*alpha**(1/2))
    u_ex['g'] = 3*c*np.cosh(scale*x)**(-2)
    u_ex['c'] *= np.exp(-1j*k*c*t)
    return u_ex['g']

# interpolate error at fixed time
def interp_err(ts, errs, t):
    try: 
        i_r = np.where(ts >= t)[0][0]
    except: 
        print("failure inside interp_err")
    if ts[i_r] == t:
        return errs[i_r]
    else:
        i_l = i_r - 1
        err = errs[i_l] + ((errs[i_r]-errs[i_l])/(ts[i_r]-ts[i_l]))*(t-ts[i_l])
        return err

# For loading IVP results (e.g., IVP under a real Fourier basis) 
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


ts = {}
errs = {}
flags = {}

for scheme in schemes:
    # load simulation data
    data_dict, flags_out = load_data(alphas, dts, scheme)
    flags[scheme] = flags_out
    errs[scheme] = {}
    ts[scheme] = {}
    for i, alpha in enumerate(alphas):
        errs[scheme][i] = {}
        ts[scheme][i] = {}
        for j, t_step in enumerate(dts):
            if flags_out[i][j]:
                print("Working on errors for scheme", scheme, i, j)
                ts_sim = data_dict[i][j]["ts_f"]
                idx_upto = np.where(ts_sim > 2*t_comp)[0][0]
                ts[scheme][i][j] = ts_sim[:idx_upto]
                errs[scheme][i][j] = np.zeros(idx_upto)
                us = data_dict[i][j]["u_f"]
                for m, t in enumerate(ts[scheme][i][j]):
                    u_ex['g'] = uex(c0, alpha, t)
                    u_sim['g'] = us[m,:]
                    u_err['g'] = u_sim['g'] - u_ex['g']
                    u_err_NL['g'] = ( -u_err*dx(u_err) ).evaluate()['g']
                    errs[scheme][i][j][m] = np.sqrt(d3.Integrate(np.abs(u_err_NL)**2, xcoord).evaluate()['g'][0])

fig13_to_plot = {}
for sch, scheme in enumerate(schemes):
    fig13_to_plot[scheme] = {}
    for i, alpha in reversed(list(enumerate(alphas))):
        dts_plot = []
        errs_NL_plot = []
        for j, t_step in enumerate(dts):
            if flags[scheme][i][j]:
                err_NL = interp_err(ts[scheme][i][j], errs[scheme][i][j], t_comp)
                errs_NL_plot.append(err_NL)
                dts_plot.append(t_step)
        fig13_to_plot[scheme][i] = {}
        fig13_to_plot[scheme][i]['xaxis'] = dts_plot
        fig13_to_plot[scheme][i]['yaxis'] = errs_NL_plot

# output processed data
np.save("processed_data/fig13-processed.npy", fig13_to_plot)
