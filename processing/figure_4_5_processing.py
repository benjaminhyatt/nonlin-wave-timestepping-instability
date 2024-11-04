"""
Processing script for Figure 4 AND Figure 5

**Assumes multiple scales analysis predictions in the L = 10 case have 
already been computed using ms_integrate_predictions.py**

(can also do processing for the L -> infty case)

Reads in IVP data from the current directory and 
multiple scales analysis predictions from the multiple scales analysis 
(stored in processed_data) and outputs fig45-processed.npy for later plotting.
"""
import numpy as np
import h5py
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)
# for integrating solvability condition with finite ansatz case
from solvability_expr_final import *
from kdv_integrate_finite_final import integrate_finite

basis = 'f' # e.g., in order to compare predictions to real Fourier IVP
ansatz = 'finite' # e.g., L = 10, the other option being 'infinite' (L -> infty)
schemes = ['sbdf1', 'sbdf2', 'rk222']

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
                data_dict[i][j]["f_flag"] = 1
            except:
                data_dict[i][j]["f_flag"] = 0

    return data_dict
    
# returns value of epsilon
def epsilon(alpha, t_step, c, L):
    FWHM = 4*float(alpha)**(1/2)*c**(-1/2)*np.log(np.sqrt(2)+1)
    T_adv = FWHM/c
    epsilon = float(t_step)/T_adv
    return epsilon

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

# compute epsilon over parameters surveyed
eps = np.zeros((alphas.shape[0], dts.shape[0]))
for i, alpha in enumerate(alphas):
    for j, t_step in enumerate(dts):
        eps[i, j] = epsilon(alpha, t_step, c0, L)

t_b_predictions = np.zeros((alphas.shape[0], dts.shape[0], len(schemes)))

# Predictions
for sch, scheme in enumerate(schemes):
    if ansatz == 'infinite':
        logger.info("Computing predictions in L -> infty limit")
        for i, alpha in enumerate(alphas):
            for j, t_step in enumerate(dts):
                t_b_pred_inf = t_b_pred_infs(c0, t_step, alpha, scheme)
                t_b_predictions[i, j, sch] = t_b_pred_inf
    
    if ansatz == 'finite':
        logger.info("Loading predictions for finite L")
    
        data_fin = np.load('processed_data/data_fin_'+scheme+'.npy', allow_pickle = True)[()]
    
        # consolidate predictions
        for i, alpha in enumerate(alphas):
            for j, t_step in enumerate(dts):
                t_pred = data_fin[i][j]['ts']
                t_b_predictions[i, j, sch] = t_pred[-1]

# IVP
fig45_to_plot = {}
for sch, scheme in enumerate(schemes):
    fig45_to_plot[scheme] = {}
    print("Loading IVP data for", scheme)
    data_dict = load_data(alphas, dts, scheme)
    for i, alpha in enumerate(alphas):
        eps_plot = []
        t_ends_a = []
        t_ends_a_pred = []
        for j, t_step in enumerate(dts):
            if basis == 'f' and data_dict[i][j]["f_flag"]:
                t_end = data_dict[i][j]["ts_f"][-1]
                eps_plot.append(eps[i, j])
                t_ends_a.append(alpha**(-1/2)*t_end)
                t_ends_a_pred.append(alpha**(-1/2)*t_b_predictions[i, j, sch])
        fig45_to_plot[scheme][i] = {}
        fig45_to_plot[scheme][i]['xaxis'] = eps_plot
        fig45_to_plot[scheme][i]['yaxis_IVP'] = t_ends_a
        fig45_to_plot[scheme][i]['yaxis_MS'] = t_ends_a_pred 

# output processed data
np.save('processed_data/fig45-processed.npy', fig45_to_plot)
