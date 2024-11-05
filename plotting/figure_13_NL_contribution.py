import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import cmasher as cmr 

##### load in processed data for plotting ##### 
fig13_to_plot = np.load('processed_data/fig13-processed.npy', allow_pickle = True)[()]

schemes = ['sbdf1', 'sbdf2', 'rk222', 'rk443']
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 

##### make 4-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600
fontsize = 8

t_mar, b_mar, l_mar, r_mar = (0.32, 0.32, 0.40, 0.12)

golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
w_pad = 0.05 * w_plot
h_pad = 0.05 * h_plot
h_cbar = h_pad

h_total = t_mar + h_pad + h_cbar + h_plot + h_pad + h_plot + b_mar
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

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot + h_pad + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad) / w_total
bottom = 1 - (t_mar + h_pad + h_cbar + h_plot + h_pad + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax4 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2, ax3, ax4]

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

for ax in [ax1, ax3]:
    ax.set_ylabel(r"$L_2(u_{\rm err} \partial_x u_{\rm err})$")

for ax in [ax1, ax2]:
    ax.set_xticklabels(())
    ax.set_ylim(1e-10, 1e3)

for ax in [ax3, ax4]:
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylim(1e-13, 1e-0)
    
for sch, ax in enumerate(axs):
    scheme = schemes[sch]
    print("maxing ax for", scheme)
    max_y = 1e-3
    for i, alpha in reversed(list(enumerate(alphas))):
        ax.scatter(fig13_to_plot[scheme][i]['xaxis'], fig13_to_plot[scheme][i]['yaxis'], s = 12, c = sm0.to_rgba(alpha), marker = "d")
        if np.max(fig13_to_plot[scheme][i]['xaxis']) > max_y:
            max_y = np.max(fig13_to_plot[scheme][i]['yaxis'])

    # trend lines
    idxlow = 0 
    if scheme == 'sbdf1':
        y = np.array(fig13_to_plot[scheme][i]['xaxis']) 
        yinit = fig13_to_plot[scheme][i]['yaxis'][0] / 1e3 
        ax.plot(y, ((yinit/y[0]**2)*y**2), linestyle = "dashed", color = "black")
    
        x_midpt_tl = np.logspace(np.log10(y[0]), np.log10(y[-1]), 3)[1]
        y_midpt_tl = np.logspace(np.log10(((yinit/y[0]**2) * y**2)[0]), np.log10(((yinit/y[0]**2) * y**2)[-1]), 3)[1]
        ax.annotate(r'$\propto \Delta t^2$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(10, -10), ha = 'center')

    elif scheme == 'sbdf2' or scheme == 'rk222':
        y = np.array(fig13_to_plot[scheme][i]['xaxis'])
        yinit = fig13_to_plot[scheme][i]['yaxis'][0] * 1e2
        idxhi = np.where((yinit/y[0]**4) * y**4 > max_y)[0][1]
        ax.plot(y[idxlow:idxhi], ((yinit/y[0]**4) * y**4)[idxlow:idxhi], linestyle = "dashed", color = "black")

        x_midpt_tl = np.logspace(np.log10((y[idxlow:idxhi])[0]), np.log10((y[idxlow:idxhi])[-1]), 3)[1]
        y_midpt_tl = np.logspace(np.log10((((yinit/y[0]**(4)) * y**(4))[idxlow:idxhi])[0]), np.log10((((yinit/y[0]**(4)) * y**(4))[idxlow:idxhi])[-1]), 3)[1]
        ax.annotate(r'$\propto \Delta t^4$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, 4), ha = 'center')

    elif scheme == 'rk443':
        y = np.array(fig13_to_plot[scheme][i]['xaxis'])
        yinit = fig13_to_plot[scheme][i]['yaxis'][0] * 1e2 
        idxhi = np.where((yinit/y[0]**6) * y**6 > max_y)[0][1]
        ax.plot(y[idxlow:idxhi], ((yinit/y[0]**6) * y**6)[idxlow:idxhi], linestyle = "dashed", color = "black")

        x_midpt_tl = np.logspace(np.log10((y[idxlow:idxhi])[0]), np.log10((y[idxlow:idxhi])[-1]), 3)[1]
        y_midpt_tl = np.logspace(np.log10((((yinit/y[0]**(6)) * y**(6))[idxlow:idxhi])[0]), np.log10((((yinit/y[0]**(6)) * y**(6))[idxlow:idxhi])[-1]), 3)[1]
        ax.annotate(r'$\propto \Delta t^6$', (x_midpt_tl, y_midpt_tl), textcoords='offset points', xytext=(-10, 4), ha = 'center')


cbar = fig.colorbar(sm, cax = cax, orientation = 'horizontal', ticks = alphas, spacing='uniform')
cax.xaxis.set_ticks_position('top')
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_fontsize(6)
cax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3f'))
cax.tick_params(which='minor', color = 'w')
cbar.ax.text(.5, 5, r'$\alpha$', ha='center', va='bottom', fontsize=8, transform=cbar.ax.transAxes)

for sch, ax in enumerate(axs):
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(5e-4, 6e-2)
    ax.text(0.03, 0.95, schemes[sch].upper(), ha = 'left', va = 'top', transform=ax.transAxes)
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())

for ax in[ax2, ax4]:
    ax.set_yticklabels(())
    ax.tick_params(axis="y", direction="in", which='both')
  
ax1.set_yticks((1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3))
ax2.set_yticks((1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3))
ax1.set_yticklabels((r'$10^{-10}$', '', '', r'$10^{-7}$', '', '', r'$10^{-4}$', '', '', r'$10^{-1}$', '', '', r'$10^2$', ''))
ax3.set_yticks((1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0))
ax4.set_yticks((1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0))
ax3.set_yticklabels((r'$10^{-13}$', '', '', r'$10^{-10}$', '', '', r'$10^{-7}$', '', '', r'$10^{-4}$', '', '', r'$10^{-1}$', ''))

plt.savefig("figure-13-NL-mag.eps")
plt.clf()
