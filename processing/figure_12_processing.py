""" 
Processing script for Figure 12

Reads in VN EVP data from current directory.
Outputs fig12-processed.npy for later plotting.
"""
import numpy as np

schemes = ["sbdf1", "sbdf2", "rk222", "rk443"]
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

# returns Re(lambda) from eigenvalue for a given timestep size
def eval_to_rate(eig, t_step):
    rate = np.log(np.abs(eig))/t_step
    return rate

fig12_to_plot = {}
for scheme in schemes:
    fig12_to_plot[scheme] = np.zeros((alphas.shape[0], dts.shape[0])) 

    # Load EVP data
    print("Reading in eigenvalues assuming only the largest eigenvalues were retained")
    evals = np.load('evals_' + scheme + '_dense_survey.npy')
       
    for i, alpha in enumerate(alphas):
        dts_plot = []
        rates_evp_plot = []
        for j, t_step in enumerate(dts):
            dts_plot.append(t_step)

            eigs = evals[i, j]
            try:
                for eig in eigs[::-1]:
                    if np.logical_not(np.isclose(eig, 0+0j)):
                        break
            except:
                eig = eigs
            rates_evp_plot.append(eval_to_rate(eig, t_step))
        
        fig12_to_plot[scheme][i] = {}
        fig12_to_plot[scheme][i]['xaxis'] = dts_plot
        fig12_to_plot[scheme][i]['yaxis'] = rates_evp_plot


# output processed data
np.save("processed_data/fig12-processed.npy", fig12_to_plot)
