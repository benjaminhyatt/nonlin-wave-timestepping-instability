"""
Plotting script for Figure 6

Reads in processed data file fig6-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import cmasher as cmr 

alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 
scheme = 'rk443'

##### load in processed data for plotting ##### 
fig6_to_plot = np.load('processed_data/fig6-processed.npy', allow_pickle = True)[()]

##### make 2-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

cmap = cmr.get_sub_cmap('viridis', 0.1, 0.85)

t_mar, b_mar, l_mar, r_mar = (0.30, 0.36, 0.32, 0.12)
golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
h_cbar = 0.05 * h_plot
h_pad = 0.05 * h_plot
w_pad = 0.1 * h_plot

h_total = t_mar + h_pad + h_cbar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + l_mar + w_plot + r_mar

fig_width = 7
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

# set up indiv ax
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2]

# set up cax
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_cbar) / h_total
width = (2. * w_plot + 1. * w_pad + 1. * l_mar) / w_total
height = h_cbar / h_total
cax = fig.add_axes([left, bottom, width, height])

norm0 = matplotlib.colors.LogNorm(vmin=alphas[0], vmax = alphas[-1])
sm0 = matplotlib.cm.ScalarMappable(norm = norm0, cmap = cmap)
dlog = 0.5 * (np.log10(alphas[1]) - np.log10(alphas[0]))
bounds = np.logspace(np.log10(alphas[0]) - dlog, np.log10(alphas[-1]) + dlog, 11)
colors = matplotlib.colors.ListedColormap([sm0.to_rgba(k) for k in alphas])
norm = matplotlib.colors.BoundaryNorm(bounds, colors.N)
sm = matplotlib.cm.ScalarMappable(norm = norm, cmap = colors)

ax = ax1
print("maxing T_decay ax")

for i, alpha in enumerate(alphas):
    ax.scatter(fig6_to_plot[i]['xaxis'], fig6_to_plot[i]['yaxis_IVP'], color = sm0.to_rgba(alpha), linewidth = 1, marker = "^")

    x_midpt_tl = 2e-2
    y_midpt_tl = np.logspace(np.log10(fig6_to_plot[i]['yaxis_MS'][0]), np.log10(fig6_to_plot[i]['yaxis_MS'][-1]), 3)[1]
    if i == 4 :
       y_midpt_tl *= 8 * 10**0
       ax.annotate(r'$\propto \epsilon^{-3}$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, -12), ha = 'center')
    ax.plot(fig6_to_plot[i]['xaxis'], fig6_to_plot[i]['yaxis_MS'], color = "black")

sim = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[4]), marker = "^", markersize = 8, linestyle = '', label = "Sim")
pred = matplotlib.lines.Line2D([], [], color = "black", linestyle = "solid", label = "MS")
ax.legend(handles=[sim, pred], loc = 'upper right', borderaxespad=0.05, fontsize = fontsize, frameon=False)

ax.text(0.05, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax.transAxes)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'$\epsilon = \frac{\Delta t \alpha^{-1/2} c_0^{3/2}}{4\ln\left(1+\sqrt{2}\right)}$')
ax.set_ylim(3e1, 5e7)
ax.set_yticks((1e2, 1e3, 1e4, 1e5, 1e6, 1e7))
ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_yticklabels(('', r'$10^3$', '', r'$10^5$', '', r'$10^7$'))
ax.set_ylabel(r'$\left(T_{\rm decay}\right)/\sqrt{\alpha}$')


ax = ax2
print("Making relative error in T_decay prediction ax")

for i, alpha in enumerate(alphas):
    y_IVP = np.array(fig6_to_plot[i]['yaxis_IVP'])[:,0]
    y_MS = np.array(fig6_to_plot[i]['yaxis_MS'])
    y_rel_err = np.abs(y_IVP - y_MS)/y_MS

    ax.scatter(np.array(fig6_to_plot[i]['xaxis']), y_rel_err, color = sm0.to_rgba(alpha), linewidth = 1, marker = "^")

idxlow = 0
y = np.array(fig6_to_plot[i]['xaxis'])
yinit = y_rel_err[0] * 8e0
try:
    idxhi = np.where((yinit/y[0]**(2)) * y**(2) > 4e0)[0][1]
except:
    idxhi = -1
ax.plot(y[idxlow:idxhi], ((yinit/y[0]**(2)) * y**(2))[idxlow:idxhi], linestyle = "dashed", color = "black")

x_midpt_tl = np.logspace(np.log10((y[idxlow:idxhi])[0]), np.log10((y[idxlow:idxhi])[-1]), 3)[1]
y_midpt_tl = np.logspace(np.log10((((yinit/y[0]**(2)) * y**(2))[idxlow:idxhi])[0]), np.log10((((yinit/y[0]**(2)) * y**(2))[idxlow:idxhi])[-1]), 3)[1]
ax.annotate(r'$\propto \epsilon^{2}$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, 4), ha = 'center')

ax.text(0.05, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax.transAxes)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'$\epsilon = \frac{\Delta t \alpha^{-1/2} c_0^{3/2}}{4\ln\left(1+\sqrt{2}\right)}$')

ax.set_ylim(1e-3, 1e1)
ax.set_yticks((1e-3, 1e-2, 1e-1, 1e0, 1e1))
ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_yticklabels((r'$10^{-3}$', '', r'$10^{-1}$', '', r'$10^{1}$'))
ax.set_ylabel(r'$\left|T_{\rm decay} - T_{\rm decay, MS}\right|\, /\, T_{\rm decay, MS}$')


cbar = fig.colorbar(sm, cax = cax, orientation = 'horizontal', ticks = alphas, spacing='uniform')
cax.xaxis.set_ticks_position('top')
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_fontsize(6)
cax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3f'))
cax.tick_params(which='minor', color = 'w')
cbar.ax.text(.5, 5, r'$\alpha$', ha='center', va='bottom', fontsize=8, transform=cbar.ax.transAxes)

plt.savefig("figure-6-t-end-rk443-2-panel.eps")
plt.clf()















