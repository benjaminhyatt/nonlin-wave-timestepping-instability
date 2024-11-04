"""
Processing script for Figure 6 (similar to that for Figures 4 and 5, but treats
the distinct case of RK443)

**Assumes multiple scales analysis predictions in the L = 10 case have 
already been computed using ms_integrate_predictions.py**

(can also do processing for the L -> infty case)

Reads in IVP data from the current directory and 
multiple scales analysis predictions from the multiple scales analysis 
(stored in processed_data) and outputs fig6-processed.npy for later plotting.
"""
import numpy as np
import h5py
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

df = 0.90 # decay fraction in L2 norm 

basis = 'f' # e.g., in order to compare predictions to real Fourier IVP
ansatz = 'finite' # e.g., L = 10, the other option being 'infinite' (L -> infty)
scheme = 'rk443'

L = 10
c0 = 0.5 
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

# load in restart flags (our RK443 simulations were long so 
# several had to be restarted from an intermediate point)
restart_flags = np.array(np.load("rk443-restart-flags.npy"))

# Load IVP data (e.g., for a real Fourier simulation)
def load_data(alphas, dts, scheme, restart_flags):
    data_dict = {}
    for i, alpha in enumerate(alphas):
        data_dict[i] = {}
        for j, t_step in enumerate(dts):
            data_dict[i][j] = {}
            if scheme == 'rk443' and not restart_flags[i, j]:
                try:
                    output_name_f = scheme + f'_output_f_{i}_{j}'
                    file_in = h5py.File(scheme + "_outputs_N_512_f/" + output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
                    data_dict[i][j]["ts_f"] = np.array(file_in['scales/sim_time'])
                    data_dict[i][j]["int_u_sq_f"] = np.array(file_in['tasks']['int_u_f_sq'])
                    data_dict[i][j]["f_flag"] = 1
                except:
                    data_dict[i][j]["f_flag"] = 0
            elif scheme == 'rk443' and restart_flags[i, j]:
                try:
                    output_name_f = scheme + f'_output_f_{i}_{j}'
                    file_in1 = h5py.File(scheme + "_outputs_N_512_f/" + output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
                    file_in2 = h5py.File(scheme + "_outputs_N_512_f/" + output_name_f + '/' + output_name_f + '_s2.h5', 'r+')
                    ts1 = np.array(file_in1['scales/sim_time'])
                    ts2 = np.array(file_in2['scales/sim_time'])
                    data_dict[i][j]["ts_f"] = np.concatenate((ts1, ts2))
                    int_u_sq_f1 = np.array(file_in1['tasks']['int_u_f_sq'])
                    int_u_sq_f2 = np.array(file_in2['tasks']['int_u_f_sq'])
                    data_dict[i][j]["int_u_sq_f"] = np.concatenate((int_u_sq_f1, int_u_sq_f2))
                    data_dict[i][j]["f_flag"] = 1 
                except:
                    data_dict[i][j]["f_flag"] = 0 
    return data_dict

# returns the (squared) KdV soliton L2 norm predicted by either ansatz 
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

# determine time at which a given threshold in decay of L2 norm is first passed, using a linear
# approximation between data points in the time series
def interp_t_l2(ts, l2s, l20, f):
    try:
        i_r = np.where(l2s < l20*f)[0][0]
        i_l = i_r - 1
        t_pred = ts[i_l] + ((ts[i_r]-ts[i_l])/(l2s[i_r]-l2s[i_l]))*(l20*f-l2s[i_l])
    except:
        logger.info("alert of interp_t exception")
        t_pred = ts[-1]
    return t_pred

# compute epsilon over parameters surveyed
eps = np.zeros((alphas.shape[0], dts.shape[0]))
for i, alpha in enumerate(alphas):
    for j, t_step in enumerate(dts):
        eps[i, j] = epsilon(alpha, t_step, c0, L)

t_b_predictions = np.zeros((alphas.shape[0], dts.shape[0]))

# Predictions
if ansatz == 'infinite':
    logger.info("Computing predictions in L -> infty limit")
    for i, alpha in enumerate(alphas):
        for j, t_step in enumerate(dts):
            t_b_pred_inf = t_b_pred_infs(c0, t_step, alpha, scheme)
            t_b_predictions[i, j] = (df**(-4) - 1)*(-1 * t_b_pred_inf) 

    logger.info("ending infinite ansatz processing tasks")

elif ansatz == 'finite':
    logger.info("Loading predictions for finite L")

    data_fin = np.load('processed_data/data_fin_'+scheme+'.npy', allow_pickle = True)[()]

    # consolidate predictions
    for i, alpha in enumerate(alphas):
        for j, t_step in enumerate(dts):
            t_pred = data_fin[i][j]['ts']
            c_pred = data_fin[i][j]['c']
            l2_pred = l2(c_pred, alpha, L, 'finite')
            l20 = l2(c0, alpha, L, ansatz)
            t_b_predictions[i, j] = interp_t_l2(t_pred, l2_pred, l20, df)

# IVP 
fig6_to_plot = {}
print("Loading IVP data for", scheme)
data_dict = load_data(alphas, dts, scheme, restart_flags)
for i, alpha in enumerate(alphas):
    eps_plot = []
    t_ends_a = []
    t_ends_a_pred = []
    for j, t_step in enumerate(dts):
        if basis == 'f' and data_dict[i][j]["f_flag"]:
            l20 = l2(c0, alpha, L, ansatz)
            t_end = interp_t_l2(data_dict[i][j]["ts_f"], np.sqrt(data_dict[i][j]["int_u_sq_f"]), l20, df)
            eps_plot.append(eps[i, j])
            t_ends_a.append(alpha**(-1/2)*t_end)
            t_ends_a_pred.append(alpha**(-1/2)*t_b_predictions[i, j])
    fig6_to_plot[i] = {}
    fig6_to_plot[i]['xaxis'] = eps_plot
    fig6_to_plot[i]['yaxis_IVP'] = t_ends_a
    fig6_to_plot[i]['yaxis_MS'] = t_ends_a_pred

# output processed data
np.save('processed_data/fig6-processed.npy', fig6_to_plot)
                                                             


