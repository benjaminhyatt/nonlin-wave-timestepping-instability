"""
This script processes eigenmodes that were found using the procedure outlined
in kdv_vn_evp.py. 

Reads in eigenvalues from EVPs computed with different N 
(expected to be saved in the current directory).

Returns the processed eigenvalues, sorted as resolved or under-resolved,
in processed_data/evals_dense_example_{i}_{j}_processed.npy 
(where i and j correspond to the particular choice of alpha and t_step) being studied).
"""
import numpy as np

# follows Chapter 7 of Boyd to accept/reject evals
# inputs should have removed inf and nans
def separate_resolved(evals_N_lo, evals_N_hi):
    
    # ensure sorted by absolute value
    idx_lo = np.abs(evals_N_lo).argsort()[::-1]
    evals_N_lo_sort = evals_N_lo[idx_lo]
    idx_hi = np.abs(evals_N_hi).argsort()[::-1]
    evals_N_hi_sort = evals_N_hi[idx_hi]

    len_lo, len_hi = (evals_N_lo_sort.shape[0], evals_N_hi_sort.shape[0])

    # calculate separations between adjacent (in abs val) modes from N_lo evals
    sigmas_lo = np.zeros(len_lo)
    sigmas_lo[0] = np.abs(evals_N_lo_sort[1] - evals_N_lo_sort[0])
    for i in range(len_lo-1):
        sigmas_lo[i] = 0.5*(np.abs(evals_N_lo_sort[i] - evals_N_lo_sort[i-1]) + np.abs(evals_N_lo_sort[i+1]-evals_N_lo_sort[i]))
    sigmas_lo[-1] = np.abs(evals_N_lo_sort[-2] - evals_N_lo_sort[-1]) 
    
    # calulate scaled differences between nearest evals from N_lo and N_hi
    idx_nearest = [np.argmin(np.abs(evals_N_lo_sort[i] - evals_N_hi_sort)/sigmas_lo[i]) for i in range(len_lo)]
    deltas = np.array([np.abs(evals_N_lo_sort[i] - evals_N_hi_sort[idx_nearest[i]])/sigmas_lo[i] for i in range(len_lo)]) 
    
    # apply threshold to inverse deltas
    drifts = 1/deltas
    drift_thresh = 1e3

    # separate 
    
    evals_resolved = evals_N_lo_sort[np.where(drifts > drift_thresh)[0]]
    evals_unresolved = evals_N_lo_sort[np.where(drifts <= drift_thresh)[0]]    
    idx_resolved = []
    idx_unresolved = []

    for eig in evals_resolved:
        idx_resolved.append(np.where(evals_N_lo_sort == eig))
    for eig in evals_unresolved:
        idx_unresolved.append(np.where(evals_N_lo_sort == eig))

    return evals_resolved, evals_unresolved, idx_resolved, idx_unresolved

# select parameters to process
alphas = np.logspace(np.log10(3e-3), np.log10(2e-2), 10) 
dts = np.logspace(np.log10(6e-4), np.log10(5e-2), 22) 
i, j = (4, 8)
alpha, t_step = (alphas[i], dts[j])

# low and high resolution evp results
N_lo, N_hi = (512, 768)

# timestepping schemes for which evp was performed
schemes = ['sbdf1', 'sbdf2', 'sbdf3', 'sbdf4', 'rk222', 'rk443']

evals_dict = {}
for scheme in schemes:

    # loading eigenvalues from evps
    evals_dict[scheme] = {}
    print('Loading N_lo results for', scheme)
    evals_dict[scheme][N_lo] = {}
    evals_dict[scheme][N_lo]['sigma'] = np.load(f'evp_out/evals_{scheme}_dense_example_{i}_{j}_{N_lo}.npy') 
    evals_dict[scheme][N_lo]['lambda'] = np.log(evals_dict[scheme][N_lo]['sigma']) / t_step
    print('Loading N_hi results for', scheme)
    evals_dict[scheme][N_hi] = {}
    evals_dict[scheme][N_hi]['sigma'] = np.load(f'evp_out/evals_{scheme}_dense_example_{i}_{j}_{N_hi}.npy')
    evals_dict[scheme][N_hi]['lambda'] = np.log(evals_dict[scheme][N_hi]['sigma']) / t_step

    # separate resolved modes
    evals_res, evals_unres, idx_res, idx_unres = separate_resolved(evals_dict[scheme][N_lo]['sigma'], evals_dict[scheme][N_hi]['sigma'])
    print('recovered', evals_res.shape[0], 'resolved modes, and', evals_unres.shape[0], 'under-resolved for scheme', scheme)

    evals_dict[scheme][N_lo]['sigma_res'] = evals_res
    evals_dict[scheme][N_lo]['sigma_unres'] = evals_unres
    evals_dict[scheme][N_lo]['lambda_res'] = evals_dict[scheme][N_lo]['lambda'][idx_res]
    evals_dict[scheme][N_lo]['lambda_unres'] = evals_dict[scheme][N_lo]['lambda'][idx_unres]

# output processed evals
np.save(f'processed_data/evals_dense_example_{i}_{j}_processed.npy', evals_dict)
