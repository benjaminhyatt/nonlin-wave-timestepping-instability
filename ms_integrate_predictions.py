"""
Script used to run numerical integrations of the nontrivial solvability conditions
for a specified IMEX scheme used to generate predictions from the multiple scales analysis. 

Outputs a time series for the effective soliton parameter c as a function of time in
processed_data/data_fin_{scheme}.npy
"""

import numpy as np
import h5py
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

# for integrating solvability condition with finite ansatz case -- note these files
# also have import dependencies of their own
from solvability_expr_final import *
from kdv_integrate_finite_final import integrate_finite

schemes = ['sbdf1', 'sbdf2', 'rk222', 'rk443']
s = 1 # (e.g., for sbdf2)
scheme = schemes[s]

L = 10
c0 = 0.5 
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

# returns value of epsilon
def epsilon(alpha, t_step, c, L): 
    FWHM = 4*float(alpha)**(1/2)*c**(-1/2)*np.log(np.sqrt(2)+1)
    T_adv = FWHM/c
    epsilon = float(t_step)/T_adv
    return epsilon

# returns hard stop time of the numerical integration
# (which is chosen to be well-beyond the expected time needed)
def T_stop(scheme, L, c0, alpha, t_step):
    if scheme == 'sbdf1':
        # blow-up time as L -> infty
        T = (35/34) * c0**(-3) * alpha * t_step**(-1)
        t_stop = 2 * T * epsilon(alpha, t_step, c0, L)
    elif scheme == 'sbdf2':
        # blow-up time as L -> infty
        T = (35/86) * c0**(-6) * alpha**2 * t_step**(-3)
        t_stop = 2 * T * (epsilon(alpha, t_step, c0, L))**3 
    elif scheme == 'rk222':
        # blow-up time as L -> infty
        T = (5005/(355045-245436*np.sqrt(2))) * c0**(-6) * alpha**2 * t_step**(-3)
        t_stop = 2 * T * (epsilon(alpha, t_step, c0, L))**3
    elif scheme == 'rk443':
        # related to the decay time as L -> infty
        T = (30030/77069) * c0**(-6) * alpha**2 * t_step**(-3)
        t_stop = 3 * T * (epsilon(alpha, t_step, c0, L))**3
    else:
        print("Not implemented", scheme)
        raise
    return t_stop


# load in expressions of nontrivial solvability conditions
if scheme == 'sbdf1':
    LHS_str, RHS_str = get_expr_sbdf1()
if scheme == 'sbdf2':
    LHS_str, RHS_str = get_expr_sbdf2()
if scheme == 'rk222':
    LHS_str, RHS_str = get_expr_rk222()
if scheme == 'rk443':
    LHS_str, RHS_str = get_expr_rk443()

# numerical integration over parameters
data_fin = {}
for i, alpha in enumerate(alphas):
    data_fin[i] = {}
    for j, t_step in enumerate(dts):
        data_fin[i][j] = {}
        logger.info("i=%i j=%i integrating solvability condition" %(i, j))
        t_stop = T_stop(scheme, L, c0, alpha, t_step)
        t_pred_fin, c_pred_fin = integrate_finite(alpha, t_step, L, c0, LHS_str, RHS_str, scheme, t_stop)

# save results
np.save('processed_data/data_fin_'+scheme+'.npy', data_fin)
