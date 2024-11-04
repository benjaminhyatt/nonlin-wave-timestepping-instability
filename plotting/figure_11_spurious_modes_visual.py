""" 
Plotting script for Figure 11

Reads in processed data file fig11-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib.pyplot as plt 

##### load in processed data for plotting ##### 
fig11_to_plot = np.load('processed_data/fig11-processed.npy', allow_pickle = True)[()]

scheme = 'rk222'

Ns = np.array((256, 384, 512)) # first dict idx
cutoffs = np.logspace(-2, -16, 15) # second dict idx

##### make 2-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

markers_M = ['1', '2', '3', 'o']
markersizes_M = [40, 40, 40, 8]
lws_M = [0.7, 0.7, 0.7, 1]
labels_M = [r'cutoff $=10^{-2}$', r'cutoff $=10^{-3}$', r'cutoff $=10^{-4}$',r'cutoff$=10^{-16}$']
colors_M = ['#dfc27d', '#bf812d', '#8c510a', '#543005']

colors_N = ['#80cdc1', '#35978f', '#01665e']

t_mar, b_mar, l_mar, r_mar = (0.12, 0.32, 0.40, 0.12)

golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
w_pad = 0.05 * w_plot
h_pad = 0.05 * h_plot

h_total = t_mar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + l_mar + w_plot + r_mar

fig_width = 7
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

##### LEFT column panels: zoomed-in scatter plots of sigma with N = 512 and 3 values of ncc cutoff #####
print("maxing ax1")
N = 512 

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1 = fig.add_axes([left, bottom, width, height])

lambdas = fig11_to_plot["ax1"]["lambda"]
cutoffs = np.array(list(lambdas[N].keys()), dtype = np.float64)
cutoffs = np.concatenate((cutoffs[cutoffs > 9e-5], cutoffs[cutoffs < 5e-16]))
for m, M in enumerate(cutoffs):
    ax1.scatter(lambdas[N][M].real, lambdas[N][M].imag, c = colors_M[m], s = markersizes_M[m], marker = markers_M[m], linewidths = lws_M[m], label = labels_M[m])

ax1.set_xlabel(r"Re($\lambda$)")
ax1.set_ylabel(r"Im($\lambda$)")
ax1.set_xlim(-1.3e0, 1.3e0)
ax1.set_ylim(-1e0, 1e0)
ax1.text(0.03, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax1.transAxes)

ax1.legend(loc = "upper right", borderaxespad=0.05, fontsize = 7, frameon=True)

##### RIGHT column panels: relative changes in sigma with respect to N and ncc cutoff #####
print("making ax2")
M_ref = cutoffs[-1]

left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2 = fig.add_axes([left, bottom, width, height])

ax2.set_xscale("log")
ax2.set_yscale("log")

for n, N in enumerate(Ns):
    ax2.plot(fig11_to_plot[n]["xaxis"], fig11_to_plot[n]["yaxis"], color = colors_N[n], marker = "o", markersize = 3, label = rf"$N = {N}$")
# emphasis
ax2.scatter(fig11_to_plot["xaxis-emph"], fig11_to_plot["yaxis-emph"], zorder = 20, color = "black", marker = "x", lw = 0.8)

ax2.legend(loc = "upper left", borderaxespad=0.05, fontsize = 7, frameon=True)

ax2.set_ylabel(r"$\|\lambda_{\rm max} - \tilde{\lambda}_{\rm max}\|_2$")
ax2.set_xlabel("NCC cutoff")
ax2.set_ylim(1e-14, 1e1)

plt.savefig("figure-11-spurious-modes.eps")
plt.clf()
