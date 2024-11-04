"""
Processing script for Figure 7

**Assumes multiple scales analysis predictions in the L = 10 case have 
already been computed using ms_integrate_predictions.py**

(can also do processing for the L -> infty case)

Reads in IVP data from the current directory and 
multiple scales analysis predictions from the multiple scales analysis 
(stored in processed_data) and outputs fig7-processed.npy for later plotting.
"""
import numpy as np
import h5py
import sys 
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

basis = 'f' # e.g., in order to compare predictions to real Fourier IVP
ansatz = 'finite' # e.g., L = 10, the other option being 'infinite' (L -> infty)
schemes = ['sbdf1', 'sbdf2', 'rk222']
frac_times = np.array((0.05, 0.5))

L = 10
c0 = 0.5 
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

# Load IVP data (e.g., for a real Fourier simulation)
def load_data(alphas, dts, scheme):
    data_dict = {}
    for i, alpha in enumerate(alphas):
        data_dict[i] = {}
        for j, t_step in enumerate(dts):
            data_dict[i][j] = {}
            try:
                output_name_f = scheme + f'_output_f_{i}_{j}'
                file_in = h5py.File(scheme + "_outputs_N_512_f/" + output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
                data_dict[i][j]["ts_f"] = np.array(file_in['scales/sim_time'])
                data_dict[i][j]["int_u_sq_f"] = np.array(file_in['tasks']['int_u_f_sq'])
                data_dict[i][j]["f_flag"] = 1
            except:
                data_dict[i][j]["f_flag"] = 0

    return data_dict

# blow up times predicted by infinite ansatz analysis
def t_b_pred_infs(c0, t_step, alpha, scheme):
    if scheme == 'sbdf1':
        return c0**(-3) * (35/34) * (alpha/t_step)
    elif scheme == 'sbdf2':
        return c0**(-6) * (35/86) * (alpha**2/t_step**3)
    elif scheme == 'rk222':
        return c0**(-6) * (5005/(355045-np.sqrt(2)*245436)) * (alpha**2/t_step**3)
    elif scheme == 'rk443':
        return -c0**(-6) * (30030/77069) * (alpha**2/t_step**3)
    else:
        print("Not a valid scheme")
        return None

# evolution of predicted by infinite ansatz analysis
def c_infs(t, c0, t_step, alpha, scheme):
    if scheme == 'sbdf1':
        return (c0**(-3) - (34/35)*(t_step/alpha)*t)**(-1/3)
    elif scheme == 'sbdf2':
        return (c0**(-6) - (86/35)*(t_step**3/alpha**2)*t)**(-1/6)
    elif scheme == 'rk222':
        return (c0**(-6) - ((355045-np.sqrt(2)*245436)/5005)*(t_step**3/alpha**2)*t)**(-1/6)
    elif scheme == 'rk443':
        return (c0**(-6) + (77069/30030)*(t_step**3/alpha**2)*t)**(-1/6)
    else:
        print("Not a valid scheme")
        return None

# returns the KdV soliton L2 norm predicted by either ansatz 
def tanh(z):
    return np.tanh(z)
def sech(z):
    return np.cosh(z)**(-1)
def l2(c, alpha, L, ansatz):
    if ansatz == "infinite":
        l2_sq = 24*c**(3/2)*alpha**(1/2)
    if ansatz == "finite":
        arg = c**(1/2) * (L/4) * alpha**(-1/2)
        arg0 = c0**(1/2) * (L/4) * alpha**(-1/2)
        l2_sq = (12/L)*(c**(3/2) * L * alpha**(1/2) * (2 + sech(arg)**2) * tanh(arg) - 12*c*alpha*tanh(arg)**2 + 12*c0*alpha*tanh(arg0)**2)
    return np.sqrt(l2_sq)

def interp_c(ts, cs, t):
    try:
        i_r = np.where(ts > t)[0][0]
        if i_r > 0:
            i_l = i_r - 1
            c_pred = cs[i_l] + ((cs[i_r]-cs[i_l])/(ts[i_r]-ts[i_l]))*(t-ts[i_l])
        else:
            c_pred = 0.5 + ((cs[i_r]-0.5)/(ts[i_r]))*(t)
    except:
        logger.info("alert of failure")
        c_pred = np.nan
    return c_pred

def interp_l2(ts, l2s, alpha, L, t, ansatz):
    try:
        i_r = np.where(ts > t)[0][0]
        if i_r > 0:
            i_l = i_r - 1
            l2_pred = l2s[i_l] + ((l2s[i_r]-l2s[i_l])/(ts[i_r]-ts[i_l]))*(t-ts[i_l])
        else:
            l2_pred = l2(0.5, alpha, L, ansatz) + ((l2s[i_r]-l2(0.5, alpha, L, ansatz))/(ts[i_r]))*(t)
    except:
        l2_pred = l2(0.5, alpha, L, ansatz)
    return l2_pred

# calculate epsilons for given parameter space
def epsilon(alpha, t_step, c, L):
    FWHM = 4*float(alpha)**(1/2)*c**(-1/2)*np.log(np.sqrt(2)+1)
    T_adv = FWHM/c
    epsilon = float(t_step)/T_adv
    return epsilon

eps = np.zeros((alphas.shape[0], dts.shape[0]))
for i, alpha in enumerate(alphas):
    for j, t_step in enumerate(dts):
        eps[i, j] = epsilon(alpha, t_step, c0, L)

t_b_predictions = np.zeros((alphas.shape[0],dts.shape[0], len(schemes)))
l2_predictions = np.zeros((alphas.shape[0], dts.shape[0], frac_times.shape[0], len(schemes)))

# Predictions
for sch, scheme in enumerate(schemes):
    if ansatz == 'infinite':
        logger.info("Computing predictions in L -> infty limit")
        for i, alpha in enumerate(alphas):
            for j, t_step in enumerate(dts):
                t_b_pred_inf = t_b_pred_infs(c0, t_step, alpha, scheme)
                t_b_predictions[i, j, sch] = t_b_pred_inf
                for k, f in enumerate(frac_times):
                    t_frac = f*t_b_predictions[i, j, sch]
                    c_predictions[i, j, k, sch] = c_infs(t_frac, c0, t_step, alpha, scheme)
                    l2_predictions[i, j, k, sch] = l2(c_infs(t_frac, c0, t_step, alpha, scheme), alpha, L, "infinite")

    
    if ansatz == 'finite':
        logger.info("Loading predictions of c for finite L (and then computing predicted L2 error)")
    
        data_fin = np.load('processed_data/data_fin_'+scheme+'.npy', allow_pickle = True)[()]
    
        # consolidate predictions
        for i, alpha in enumerate(alphas):
            for j, t_step in enumerate(dts):
                t_pred = data_fin[i][j]['ts']
                c_pred = data_fin[i][j]['c']
                l2_pred = l2(c_pred, alpha, L, 'finite')

                t_b_predictions[i, j, sch] = t_pred[-1]
                for k, f in enumerate(frac_times):
                    t_frac = f*t_b_predictions[i, j, sch]
                    l2_predictions[i, j, k, sch] = interp_l2(t_pred, l2_pred, alpha, L, t_frac, "finite")

# IVP 
fig7_to_plot = {}

for sch, scheme in enumerate(schemes):
    fig7_to_plot[scheme] = {}

    print("Loading IVP data for", scheme)
    data_dict = load_data(alphas, dts, scheme)

    idxlow = 0 
    max_y = 1e-3

    for k, f in enumerate(frac_times):
        fig7_to_plot[scheme][k] = {}
         
        for i, alpha in enumerate(alphas):
            eps_plot = []
            l2_errors = []
            for j, t_step in enumerate(dts):
                if basis == 'f' and data_dict[i][j]["f_flag"]:
                    eps_plot.append(eps[i, j])
                    t_frac = f * data_dict[i][j]["ts_f"][-1]
                    l2_f = interp_l2(data_dict[i][j]["ts_f"], np.sqrt(data_dict[i][j]["int_u_sq_f"]), alpha, L, t_frac, ansatz)        
                    l2_f_pred = l2_predictions[i, j, k, sch]
                    l2_err = np.abs(l2_f - l2_f_pred)/l2_f_pred
                    l2_errors.append(l2_err)
                    if l2_err > max_y:
                        max_y = l2_err            
            fig7_to_plot[scheme][k][i] = {}
            fig7_to_plot[scheme][k][i]['xaxis'] = eps_plot
            fig7_to_plot[scheme][k][i]['yaxis'] = l2_errors

    fig7_to_plot[scheme]['max_y'] = max_y

# output processed data
np.save('processed_data/fig7-processed.npy', fig7_to_plot)








