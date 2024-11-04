"""
Processing script for Figure 10

Reads in processed VN EVP data from inside of processed_data that has been 
processed by process_evals.py to sort-out/reject under-resolved modes. 

Outputs fig10-processed.npy for later plotting.

Options: 
    i (Line 14) index array of dispersion parameter values (alpha)
    j (Line 14) index into array of timestep sizes
"""
import numpy as np

# select parameters
N = 512 
i, j = (4, 8)
# timestepping schemes
schemes = ['sbdf1', 'sbdf2', 'sbdf3', 'sbdf4', 'rk222', 'rk443']

# load processed eigenvalues
evals_dict = np.load(f'processed_data/evals_dense_example_{i}_{j}_processed.npy', allow_pickle = True)[()]

fig10_to_plot = {}
for s, scheme in enumerate(schemes):
    fig10_to_plot[s] = {}

    sigmas = evals_dict[scheme][N]['sigma_res']
    sigmas_u = sigmas[np.abs(sigmas) > 1.]
    sigmas_s = sigmas[np.abs(sigmas) <= 1.]
 
    # look for fastest growing mode
    maxidx_s = np.where(np.abs(sigmas_u) == np.max(np.abs(sigmas)))
    fig10_to_plot[s]['maxidx_s'] = maxidx_s
    # check for complex conjugate
    try:
        maxidx_s_cc = np.where(np.isclose(sigmas_u.imag, -1 * sigmas_u[maxidx_s].imag, atol = 1e-14))    
        if sigmas_u[maxidx_s_cc]: 
            print('identified sigma cc pair of fastest growing modes', sigmas_u[maxidx_s], sigmas_u[maxidx_s_cc])
            fig10_to_plot[s]['maxidx_s_cc'] = maxidx_s_cc
        else:
            print('did not find a sigma cc of fastest growing mode')
            fig10_to_plot[s]['maxidx_s_cc'] = False
    except:
        print('did not find a sigma cc of fastest growing mode')

    # (resolved stable modes)
    fig10_to_plot[s]['sigmas_s_re'] = sigmas_s.real
    fig10_to_plot[s]['sigmas_s_im'] = sigmas_s.imag
    # (resolved unstable modes)
    fig10_to_plot[s]['sigmas_u_re'] = sigmas_u.real
    fig10_to_plot[s]['sigmas_u_im'] = sigmas_u.imag


    lambdas = evals_dict[scheme][N]['lambda_res']
    lambdas_u = lambdas[lambdas.real > 0.]
    lambdas_s = lambdas[lambdas.real <= 0.]
    maxidx_l = np.where(lambdas_u.real == np.max(lambdas_u.real))
    fig10_to_plot[s]['maxidx_l'] = maxidx_l
    # check for complex conjugate
    try:
        maxidx_l_cc = np.where(np.isclose(lambdas_u.imag, -1 * lambdas_u[maxidx_l].imag, atol = 1e-14))
        if lambdas_u[maxidx_s_cc]:
            print('identified lambda cc pair of fastest growing modes', lambdas_u[maxidx_l], lambdas_u[maxidx_l_cc])
            fig10_to_plot[s]['maxidx_l_cc'] = maxidx_l_cc
        else:
            print('did not find a lambda cc of fastest growing mode')
            fig10_to_plot[s]['maxidx_l_cc'] = False
    except:
        print('did not find a lambda cc of fastest growing mode')

    # (resolved stable modes)
    fig10_to_plot[s]['lambdas_s_re'] = lambdas_s.real
    fig10_to_plot[s]['lambdas_s_im'] = lambdas_s.imag
    # (resolved unstable modes)
    fig10_to_plot[s]['lambdas_u_re'] = lambdas_u.real
    fig10_to_plot[s]['lambdas_u_im'] = lambdas_u.imag

# output processed data
np.save('processed_data/fig10-processed.npy', fig10_to_plot)

