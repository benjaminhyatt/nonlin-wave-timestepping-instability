""" 
Processing script for Figure 11

Reads in VN EVP data from current directory (did not apply process_evals.py
to this dataset).

Outputs fig11-processed.npy for later plotting.

Options: 
    scheme (Line 19) which scheme's spectra to be examined (currently RK222)
    i (Line 25) index array of dispersion parameter values (alpha)
    j (Line 25) index into array of timestep sizes
"""
import numpy as np
import matplotlib.pyplot as plt 

scheme = 'rk222'

c = 0.5
L = 10
alphas = np.logspace(np.log10(3e-3), np.log10(2e-2), 10)
dts = np.logspace(np.log10(6e-4), np.log10(5e-2), 22)
i, j = (4, 8)
alpha, t_step = (alphas[i], dts[j])

Ns = np.array((256, 384, 512)) # first dict idx
cutoffs = np.logspace(-2, -16, 15) # second dict idx

fig11_to_plot = {}

print("ax1 data")

fig11_to_plot["ax1"] = {}
fig11_to_plot["ax1"]["sigma"] = np.load("evals_" + scheme + "_dense_example_4_8.npy", allow_pickle = True)[()]
fig11_to_plot["ax1"]["lambda"] = {}
for n, N in enumerate(Ns):
    fig11_to_plot["ax1"]["lambda"][N] = {}
    cutoffs = np.array(list(fig11_to_plot["ax1"]["sigma"][N].keys()), dtype = np.float64)
    for m, M in enumerate(cutoffs):
        fig11_to_plot["ax1"]["lambda"][N][M] = np.log(fig11_to_plot["ax1"]["sigma"][N][M]) / t_step


print("ax2 data")
lambdas = fig11_to_plot["ax1"]["lambda"]
M_ref = cutoffs[-1]

for n, N in enumerate(Ns):
    euc_err = []
    cutoffs = list(lambdas[N].keys())
    for m, M in enumerate(cutoffs[:-1]):
        lambdas_u = lambdas[N][M][lambdas[N][M].real > 0.] 
        lambdas_u_ref = lambdas[N][M_ref][lambdas[N][M_ref].real > 0.] 
        lambda_m = lambdas_u[lambdas_u.real == np.max(lambdas_u.real)]
        if lambda_m.shape[0] > 1:
            lambda_m = lambda_m[0]
        lambda_ref = lambdas_u_ref[lambdas_u_ref.real == np.max(lambdas_u_ref.real)]
        re_err = np.abs(lambda_m.real - lambda_ref.real)
        im_err_1 = np.abs(lambda_m.imag - lambda_ref.imag)
        im_err_2 = np.abs(-lambda_m.imag - lambda_ref.imag)
        im_err = np.min((im_err_1, im_err_2))
        euc = np.sqrt(re_err**2 + im_err**2)
        euc_err.append(euc)

    fig11_to_plot[n] = {}
    fig11_to_plot[n]["xaxis"] = cutoffs[:-1]
    fig11_to_plot[n]["yaxis"] = euc_err

fig11_to_plot["xaxis-emph"] = cutoffs[:3]
fig11_to_plot["yaxis-emph"] = euc_err[:3]

# output processed data
np.save('processed_data/fig11-processed.npy', fig11_to_plot)
