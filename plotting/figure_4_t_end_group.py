"""
Plotting script for Figure 4

Reads in processed data file fig45-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import cmasher as cmr 

schemes = ['sbdf1', 'sbdf2', 'rk222']
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

##### load in processed data for plotting ##### 
fig45_to_plot = np.load('processed_data/fig45-processed.npy', allow_pickle = True)[()]

##### make 3-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

cmap = cmr.get_sub_cmap('viridis', 0.1, 0.85)

t_mar, b_mar, l_mar, r_mar = (0.30, 0.44, 0.44, 0.12)
golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
h_cbar = 0.05 * h_plot
h_pad = 0.05 * h_plot
w_pad = 0.1 * h_plot

h_total = t_mar + h_pad + h_cbar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + w_plot + w_pad + w_plot + r_mar

fig_width = 7
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

# set up indiv ax
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

left = (l_mar + 2. * w_plot + 2. * w_pad) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2, ax3]

# set up cax
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_cbar) / h_total
width = (3. * w_plot + 2. * w_pad) / w_total
height = h_cbar / h_total 
cax = fig.add_axes([left, bottom, width, height])

norm0 = matplotlib.colors.LogNorm(vmin=alphas[0], vmax = alphas[-1])
sm0 = matplotlib.cm.ScalarMappable(norm = norm0, cmap = cmap)
dlog = 0.5 * (np.log10(alphas[1]) - np.log10(alphas[0]))
bounds = np.logspace(np.log10(alphas[0]) - dlog, np.log10(alphas[-1]) + dlog, 11)
colors = matplotlib.colors.ListedColormap([sm0.to_rgba(k) for k in alphas])
norm = matplotlib.colors.BoundaryNorm(bounds, colors.N)
sm = matplotlib.cm.ScalarMappable(norm = norm, cmap = colors)

for sch, ax in enumerate(axs):
    scheme = schemes[sch]

    print("maxing ax for", scheme)
    for i, alpha in enumerate(alphas):
        ax.scatter(fig45_to_plot[scheme][i]['xaxis'], fig45_to_plot[scheme][i]['yaxis_IVP'], color = sm0.to_rgba(alpha), linewidth = 1, marker = "^")
        
        x_midpt_tl = fig45_to_plot[scheme][i]['xaxis'][13]
        y_midpt_tl = np.logspace(np.log10(fig45_to_plot[scheme][i]['yaxis_MS'][0]), np.log10(fig45_to_plot[scheme][i]['yaxis_MS'][-1]), 3)[1] 
        if i == 4 and scheme == 'sbdf1':
            y_midpt_tl_0 = y_midpt_tl * 5 * 10**1
            ax.annotate(r'$\propto \epsilon^{-1}$', (x_midpt_tl, y_midpt_tl_0), textcoords='offset points', xytext=(-10, -12), ha = 'center') 
        if i == 4 and (scheme == 'sbdf2' or scheme == 'rk222'):
            y_midpt_tl_0 = y_midpt_tl * 2 * 10**1
            ax.annotate(r'$\propto \epsilon^{-3}$', (x_midpt_tl, y_midpt_tl_0), textcoords='offset points', xytext=(-10, -12), ha = 'center')
        
        ax.plot(fig45_to_plot[scheme][i]['xaxis'], fig45_to_plot[scheme][i]['yaxis_MS'], color = "black")
        
cbar = fig.colorbar(sm, cax = cax, orientation = 'horizontal', ticks = alphas, spacing='uniform') 
cax.xaxis.set_ticks_position('top')
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_fontsize(6)
cax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3f'))
cax.tick_params(which='minor', color = 'w')
cbar.ax.text(.5, 5, r'$\alpha$', ha='center', va='bottom', fontsize=8, transform=cbar.ax.transAxes)

sim = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[4]), marker = "^", markersize = 8, linestyle = '', label = "Sim")
pred = matplotlib.lines.Line2D([], [], color = "black", linestyle = "solid", label = "MS")
ax1.legend(handles=[sim, pred], loc = 'upper right', borderaxespad=0.05, fontsize = fontsize, frameon=False)


for sch, ax in enumerate(axs):
    scheme = schemes[sch]
    ax.text(0.05, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax.transAxes)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r'$\epsilon = \frac{\Delta t \alpha^{-1/2} c_0^{3/2}}{4\ln\left(1+\sqrt{2}\right)}$')
    ax.set_ylim(5e0, 5e8)
    ax.set_yticks((1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    
# non-left plots
for ax in [ax2, ax3]:
    ax.set_yticklabels(())
    ax.tick_params(axis="y", direction="in", which='both')

ax1.set_yticklabels(('', r'$10^{2}$', '', r'$10^{4}$', '', r'$10^{6}$', '', r'$10^{8}$'))
ax1.set_ylabel(r'$\left(T_{\rm blowup}\right)/\sqrt{\alpha}$')

plt.savefig("figure-4-t-ends-group.eps")
plt.clf()
