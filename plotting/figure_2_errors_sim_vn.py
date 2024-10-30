"""
Plotting script for Figure 2

Reads in processed data file fig2-processed.npy
from processed_data folder (which should be within current working directory)

"""

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import publication_settings
import cmasher as cmr 

schemes = ["sbdf1", "sbdf2", "rk222", "rk443", "sbdf3", "sbdf4"]

##### load in processed data for plotting ##### 
fig2_to_plot = np.load('processed_data/fig2-processed.npy', allow_pickle = True)[()]
i, j = fig2_to_plot['(i, j)']


##### make 6-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8
plt.rcParams['figure.dpi'] = 600
fontsize = 8

t_mar, b_mar, l_mar, r_mar = (0.12, 0.36, 0.44, 0.12)

golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
w_pad = 0.05 * w_plot
h_pad = w_pad

h_total = t_mar + h_plot + h_pad + h_plot + b_mar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + w_plot + r_mar

fig_width = 5.5
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

##### construct axs #####

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2 = fig.add_axes([left, bottom, width, height])

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + h_pad + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad) / w_total
bottom = 1 - (t_mar + h_plot + h_pad + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax4 = fig.add_axes([left, bottom, width, height])

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + h_pad + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax5 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad) / w_total
bottom = 1 - (t_mar + h_plot + h_pad + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax6 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2, ax3, ax4, ax5, ax6]

vncol = matplotlib.colors.to_rgba('#d95f02')
simcol = matplotlib.colors.to_rgba('#7570b3')

# plot axes
for sch, ax in enumerate(axs):
    scheme = schemes[sch]
    print("making ax for", scheme)

    ax.plot(fig2_to_plot[scheme]['ivp_xaxis'], fig2_to_plot[scheme]['ivp_yaxis'], linewidth = 2, color = simcol)
    ax.plot(fig2_to_plot[scheme]['evp_xaxis'], fig2_to_plot[scheme]['evp_yaxis'], linewidth = 2, color = vncol, linestyle = "--")

actual = matplotlib.lines.Line2D([], [], color = simcol, linewidth = 2, label = "Sim")
vn = matplotlib.lines.Line2D([], [], color = vncol, linewidth = 2, linestyle = "--", label = "VN")
ax1.legend(handles=[actual, vn], loc = (0.65, 0.725), fontsize = fontsize, frameon=False)

for sch, ax in enumerate(axs):
    ax.set_yscale("log")
    ax.text(0.03, 0.95, schemes[sch].upper(), ha = 'left', va = 'top', transform=ax.transAxes)
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator((0, 0.25, 0.5, 0.75, 1.0)))
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(n=2))

for ax in [ax1, ax3, ax5]:
    ax.set_ylabel(r"$L_2(u_{\rm sim}-u_{\rm ex})$")

for ax in [ax3, ax4, ax5, ax6]:
    ax.set_xticklabels(('0.0', '', r'0.5$T_{\rm end}$', '', r'$T_{\rm end}$' ))
    ax.set_xlabel(r'$t$')
    
for ax in [ax2, ax4, ax6]:
    ax.set_yticklabels(())
    ax.tick_params(axis="y", direction="in", which='both')

for ax in [ax1, ax2]:
    ax.set_xticklabels(())
    ax.set_yticks((1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2))
    ax.set_ylim(bottom = 1e-4, top = 1e2)

for ax in [ax3, ax4]:
    ax.set_yticks((1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2))
    ax.set_ylim(bottom = 1e-8, top = 1e2)

for ax in [ax5, ax6]:
    ax.set_ylim(bottom = 1e-4, top = 1e2)
    ax.set_yticks((1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2))

ax1.set_yticklabels((r'$10^{-4}$', '', r'$10^{-2}$', '', r'$10^{0}$', '', r'$10^2$'))
ax3.set_yticklabels((r'$10^{-8}$', '', r'$10^{-6}$', '', r'$10^{-4}$', '', r'$10^{-2}$', '', r'$10^{0}$', '', r'$10^2$'))
ax5.set_yticklabels((r'$10^{-4}$', '', r'$10^{-2}$', '', r'$10^{0}$', '', r'$10^2$'))

xlims5 = ax5.get_xlim()
ax5.set_xlim(xlims5)

# exponential rate labels
rotations = []
for sch, ax in enumerate(axs):
    scheme = schemes[sch]

    lamb = fig2_to_plot[scheme]['rate']
    ys_evp = fig2_to_plot[scheme]['evp_yaxis']   
 
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    if ax == ax5:
        x2 = 0.22 * xlims[1]
    elif ax == ax6:
        x2 = 0.57 * xlims[1]
    else:
        x2 = 1.
    x1 = 0.
    dx_ax = (x2 - x1) / (xlims[1] - xlims[0])  
    if ax == ax5:
        y2 = ylims[1]
    elif ax == ax6:
        y2 = ylims[1]
    else:
        y2 = np.max(ys_evp)
    y1 = np.min(ys_evp)
    dy_ax = (np.log10(y2) - np.log10(y1)) / (np.log10(ylims[1]) - np.log10(ylims[0]))
    dx_ax /= golden_mean
    angle_rad = np.arctan(dy_ax / dx_ax)
    angle_deg = (180/np.pi) * angle_rad
    rotations.append(angle_deg)

ax1.text(0.375, 7e-4, r"$\propto e^{%s t}$"%format(fig2_to_plot['sbdf1']['rate'], '.4f'), color = vncol, rotation = rotations[0])
ax2.text(0.375, 5e-4, r"$\propto e^{%s t}$"%format(fig2_to_plot['sbdf2']['rate'], '.7f'), color = vncol, rotation = rotations[1])
ax3.text(0.375, 6e-6, r"$\propto e^{%s t}$"%format(fig2_to_plot['rk222']['rate'], '.7f'), color = vncol, rotation = rotations[2])
ax4.text(0.375, 7e-7, r"$\propto e^{%s t}$"%format(fig2_to_plot['rk443']['rate'], '.8f'), color = vncol, rotation = rotations[3])
ax5.text(-0.01, 1e-2, r"$\propto e^{%s t}$"%format(fig2_to_plot['sbdf3']['rate'], '.3f'), color = vncol, rotation = rotations[4])
ax6.text(0.19, 8e-2, r"$\propto e^{%s t}$"%format(fig2_to_plot['sbdf4']['rate'], '.3f'), color = vncol, rotation = rotations[5])

plt.savefig(f"figure-2-evp-ivp-comparison_{i}_{j}.eps")
plt.clf()
