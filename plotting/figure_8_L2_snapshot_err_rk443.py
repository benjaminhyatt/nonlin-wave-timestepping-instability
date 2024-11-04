"""
Plotting script for Figure 8

Reads in processed data file fig8-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import cmasher as cmr 

frac_times = np.array((0.05, 0.5))
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 
scheme = 'rk443'

##### load in processed data for plotting ##### 
fig8_to_plot = np.load('processed_data/fig8-processed.npy', allow_pickle = True)[()]

##### make figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

cmap = cmr.get_sub_cmap('viridis', 0.1, 0.85)

t_mar, b_mar, l_mar, r_mar = (0.30, 0.36, 0.36, 0.12)
golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
h_cbar, w_cbar = (0.05 * h_plot, 2. * w_plot)
h_pad = 0.05 * h_plot

h_total = t_mar + h_pad + h_cbar + h_plot + b_mar
w_total = l_mar + w_plot + r_mar

fig_width = 3.5
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

# set up ax
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax = fig.add_axes([left, bottom, width, height])

# set up cax
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_cbar) / h_total
width = (w_plot) / w_total
height = h_cbar / h_total
cax = fig.add_axes([left, bottom, width, height])

norm0 = matplotlib.colors.LogNorm(vmin=alphas[0], vmax = alphas[-1])
sm0 = matplotlib.cm.ScalarMappable(norm = norm0, cmap = cmap)
dlog = 0.5 * (np.log10(alphas[1]) - np.log10(alphas[0]))
bounds = np.logspace(np.log10(alphas[0]) - dlog, np.log10(alphas[-1]) + dlog, 11)
colors = matplotlib.colors.ListedColormap([sm0.to_rgba(k) for k in alphas])
norm = matplotlib.colors.BoundaryNorm(bounds, colors.N)
sm = matplotlib.cm.ScalarMappable(norm = norm, cmap = colors)

for k, f in enumerate(frac_times):
    for i, alpha in enumerate(alphas):
        if k == 0:
            ax.scatter(fig8_to_plot[k][i]['xaxis'], fig8_to_plot[k][i]['yaxis'], color = sm0.to_rgba(alpha), linewidth = 1, marker = "1")
        if k == 1:
            ax.scatter(fig8_to_plot[k][i]['xaxis'], fig8_to_plot[k][i]['yaxis'], color = sm0.to_rgba(alpha), linewidth = 1, marker = "2")

idxlow = 0
y = np.array(fig8_to_plot[k][i]['xaxis'])
yinit = fig8_to_plot[k][i]['yaxis'][0] * 4e0 
try:
    idxhi = np.where((yinit/y[0]**2) * y**2 > 6e-3)[0][1]
except:
    idxhi = -1
ax.plot(y[idxlow:idxhi], ((yinit/y[0]**2) * y**2)[idxlow:idxhi], linestyle = "dashed", color = "black")

x_midpt_tl = np.logspace(np.log10((y[idxlow:idxhi])[0]), np.log10((y[idxlow:idxhi])[-1]), 3)[1]
y_midpt_tl = np.logspace(np.log10((((yinit/y[0]**(2)) * y**(2))[idxlow:idxhi])[0]), np.log10((((yinit/y[0]**(2)) * y**(2))[idxlow:idxhi])[-1]), 3)[1]
ax.annotate(r'$\propto \epsilon^2$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, 4), ha = 'center')


cbar = fig.colorbar(sm, cax = cax, orientation = 'horizontal', ticks = alphas, spacing='uniform')
cax.xaxis.set_ticks_position('top')
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_fontsize(6)
cax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3f'))
cax.tick_params(which='minor', color = 'w')
cbar.ax.text(.5, 5, r'$\alpha$', ha='center', va='bottom', fontsize=8, transform=cbar.ax.transAxes)

t0label = r'$t/T_{\rm decay} =$' + str(frac_times[0])
t1label = r'$t/T_{\rm decay} =$' + str(frac_times[1])
t0 = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[4]), marker = "1", markersize = 10, linestyle = '', label = t0label)
t1 = matplotlib.lines.Line2D([], [], color = sm0.to_rgba(alphas[4]), marker = "2", markersize = 10, linestyle = '', label = t1label)
ax.legend(handles=[t1, t0], loc = 'lower right', borderaxespad=0.05, fontsize = fontsize, frameon=False)

ax.text(0.05, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax.transAxes)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'$\epsilon = \frac{\Delta t \alpha^{-1/2} c_0^{3/2}}{4\ln\left(1+\sqrt{2}\right)}$')
ax.set_ylim(1e-5, 1e-2)
ax.set_yticks((1e-5, 1e-4, 1e-3, 1e-2))
ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())

ax.set_yticklabels((r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$'))
ax.set_ylabel(r'$\left|L_2(u) - L_2(u_{\rm MS})\right|\, /\, L_2(u_{\rm MS})$')

plt.savefig("figure-8-L2-rel-err-wrt-time-rk443.eps")
plt.clf()


