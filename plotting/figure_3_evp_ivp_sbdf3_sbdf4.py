"""
Plotting script for Figure 3

Reads in processed data file fig3-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import cmasher as cmr 

schemes = ["sbdf3", "sbdf4"]
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10)
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22)

##### load in processed data for plotting ##### 
fig3_to_plot = np.load('processed_data/fig3-processed.npy', allow_pickle = True)[()]

##### make 2-panel figure ##### 
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

t_mar, b_mar, l_mar, r_mar = (0.30, 0.30, 0.36, 0.12)

golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
w_pad = 0.05 * w_plot
h_pad = 0.05 * h_plot
h_cbar = h_pad

h_total = t_mar + h_pad + h_cbar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + w_plot + r_mar

fig_width = 5.5
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

##### construct axs #####

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2]

##### construct cax #####
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_cbar) / h_total
width = (2. * w_plot + w_pad) / w_total
height = h_cbar / h_total
cax = fig.add_axes([left, bottom, width, height])

cmap = cmr.get_sub_cmap('viridis', 0.1, 0.85)
norm0 = matplotlib.colors.LogNorm(vmin=alphas[0], vmax = alphas[-1])
sm0 = matplotlib.cm.ScalarMappable(norm = norm0, cmap = cmap)
dlog = 0.5 * (np.log10(alphas[1]) - np.log10(alphas[0]))
bounds = np.logspace(np.log10(alphas[0]) - dlog, np.log10(alphas[-1]) + dlog, 11)
colors = matplotlib.colors.ListedColormap([sm0.to_rgba(k) for k in alphas])
norm = matplotlib.colors.BoundaryNorm(bounds, colors.N)
sm = matplotlib.cm.ScalarMappable(norm = norm, cmap = colors)

for sch, ax in enumerate(axs):
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta t$")
    
    ax.set_yscale("log")

    ax.text(0.03, 0.95, schemes[sch].upper(), ha = 'left', va = 'top', transform=ax.transAxes)

for sch, ax in enumerate(axs):
    scheme = schemes[sch]
    print("maxing ax for", scheme)
    for i, alpha in enumerate(alphas):
        # IVP
        ax.scatter(fig3_to_plot[scheme]['IVP' + str(i)]['IVP_xaxis'], fig3_to_plot[scheme]['IVP' + str(i)]['IVP_yaxis'], s = 32, c = sm0.to_rgba(alpha), marker = "x", linewidths = 1)
    for i, alpha in reversed(list(enumerate(alphas))):
        # EVP
        ax.scatter(fig3_to_plot[scheme]['EVP' + str(i)]['EVP_xaxis'], fig3_to_plot[scheme]['EVP' + str(i)]['EVP_yaxis'], s = 8, c = sm0.to_rgba(alpha), marker = "o")
    # Trend line
    y = fig3_to_plot[scheme]['trend_y']
    yinit, idxlow, idxhi, x_midpt_tl, y_midpt_tl = fig3_to_plot[scheme]['trend_misc']
    ax.plot(y[idxlow:idxhi], ((yinit/y[0]**(-1)) * y**(-1))[idxlow:idxhi], linestyle = "dashed", color = "black")
    ax.annotate(r'$\propto \Delta t ^{-1}$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, -12), ha = 'center')

# additional cax and ax formatting
cbar = fig.colorbar(sm, cax = cax, orientation = 'horizontal', ticks = alphas, spacing='uniform') 
cax.xaxis.set_ticks_position('top')
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_fontsize(6)
cax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3f'))
cax.tick_params(which='minor', color = 'w')
cbar.ax.text(.5, 5, r'$\alpha$', ha='center', va='bottom', fontsize=8, transform=cbar.ax.transAxes)

cross = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[7]), marker = "x", markersize = 6, mew = 1.0, linestyle = '', label = 'Sim')
circ = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[0]), marker = "o", markersize = 6, linestyle = '', label = 'VN')
ax1.legend(handles=[cross, circ], loc = 'upper right', borderaxespad=0.05, fontsize = fontsize, frameon=False)

ax1.set_ylabel(r"Re$\left(\lambda\right)$")

for sch, ax in enumerate(axs):
    ylims = (2.5e-1, 1.5e3)
    ax.set_ylim(ylims)
    ax.set_yticks((1e0, 1e1, 1e2, 1e3))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax2.set_yticklabels(())
ax2.tick_params(axis="y", direction="in", which='both')

plt.savefig("figure-3-evp-ivp-survey-rate-comparison.eps")
plt.clf()
