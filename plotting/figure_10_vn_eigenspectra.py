"""
Plotting script for Figure 10

Reads in processed data file fig10-processed.npy
from processed_data folder (which should be within current working directory)
"""
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

schemes = ['sbdf1', 'sbdf2', 'sbdf3', 'sbdf4', 'rk222', 'rk443']

##### load in processed data for plotting ##### 
fig10_to_plot = np.load('processed_data/fig10-processed.npy', allow_pickle = True)[()]

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

t_mar, b_mar, l_mar, r_mar = (0.08, 0.36, 0.44, 0.18)

h_plot, w_plot = (1., 1.)
w_pad = 0.1 * w_plot

h_total = t_mar + h_plot + b_mar + h_plot + b_mar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar + w_plot + r_mar


fig_width = 7
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

# sbdf1 sigma
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1s = fig.add_axes([left, bottom, width, height])

# sbdf1 lambda
left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1l = fig.add_axes([left, bottom, width, height])

# sbdf2 sigma
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2s = fig.add_axes([left, bottom, width, height])

# sbdf2 lambda
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2l = fig.add_axes([left, bottom, width, height])

# sbdf3 sigma
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3s = fig.add_axes([left, bottom, width, height])

# sbdf3 lambda
left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3l = fig.add_axes([left, bottom, width, height])

# sbdf4 sigma
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax4s = fig.add_axes([left, bottom, width, height])

# sbdf4 lambda
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax4l = fig.add_axes([left, bottom, width, height])

# rk222 sigma
left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax5s = fig.add_axes([left, bottom, width, height])

# rk222 lambda
left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax5l = fig.add_axes([left, bottom, width, height])

# rk443 sigma
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax6s = fig.add_axes([left, bottom, width, height])

# rk443 lambda
left = (l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax6l = fig.add_axes([left, bottom, width, height])

axs = [ax1s, ax1l, ax2s, ax2l, ax3s, ax3l, ax4s, ax4l, ax5s, ax5l, ax6s, ax6l]

for s, scheme in enumerate(schemes):
    print("plotting axes for", scheme)
    ax1 = axs[int(2*s)]
    ax2 = axs[int(2*s)+1]

    # sigma plot
    
    # draw unit circle
    ax1.add_patch(Circle((0., 0.), 1., linestyle='--', fill = False, color = 'black'))

    # scatter resolved eigenvalues
    ax1.scatter(fig10_to_plot[s]['sigmas_s_re'], fig10_to_plot[s]['sigmas_s_im'], label = r'$|\sigma| < 1$', s = 1, c = '#1b9e77')
    ax1.scatter(fig10_to_plot[s]['sigmas_u_re'], fig10_to_plot[s]['sigmas_u_im'], label = r'$|\sigma| > 1$', s = 10, c = '#7570b3')
    # highlight fastest growing mode(s)
    try:
        maxidx_s = fig10_to_plot[s]['maxidx_s']
        maxidx_s_cc = fig10_to_plot[s]['maxidx_s_cc']
        ax1.scatter([fig10_to_plot[s]['sigmas_u_re'][maxidx_s], fig10_to_plot[s]['sigmas_u_re'][maxidx_s_cc]], [fig10_to_plot[s]['sigmas_u_im'][maxidx_s], fig10_to_plot[s]['sigmas_u_im'][maxidx_s_cc]], label = r'$|\sigma|_{\rm max}$', s = 30, c = '#d95f02', marker = 'x', linewidth = 1.)
    except:
        maxidx_s = fig10_to_plot[s]['maxidx_s']
        ax1.scatter(fig10_to_plot[s]['sigmas_u_re'][maxidx_s], fig10_to_plot[s]['sigmas_u_im'][maxidx_s], label = r'$|\sigma|_{\rm max}$', s = 30, c = '#d95f02', marker = 'x', linewidth = 1.)

    ax1.set_xlabel(r'Re($\sigma$)')
    ax1.set_ylabel(r'Im($\sigma$)')

    ax1.set_xlim(-1.35, 1.35)
    ax1.set_ylim(-1.35, 1.35)
    if scheme == 'sbdf1':
        ax1.legend(loc = 'lower left', borderaxespad=0.05, fontsize = 6, frameon=True)
    ax1.text(0.03, 0.95, scheme.upper(), ha = 'left', va = 'top', transform=ax1.transAxes)

    # lambda plot
    
    # scatter resolved eigenvalues
    ax2.scatter(fig10_to_plot[s]['lambdas_s_re'], fig10_to_plot[s]['lambdas_s_im'], label = r'Re$(\lambda) < 0$', s = 1, c = '#1b9e77')
    ax2.scatter(fig10_to_plot[s]['lambdas_u_re'], fig10_to_plot[s]['lambdas_u_im'], label = r'Re$(\lambda) > 0$', s = 10, c = '#7570b3')
    # highlight fastest growing mode(s)
    try:
        maxidx_l = fig10_to_plot[s]['maxidx_l']
        maxidx_l_cc = fig10_to_plot[s]['maxidx_l_cc']
        ax2.scatter([fig10_to_plot[s]['lambdas_u_re'][maxidx_l], fig10_to_plot[s]['lambdas_u_re'][maxidx_l_cc]], [fig10_to_plot[s]['lambdas_u_im'][maxidx_l], fig10_to_plot[s]['lambdas_u_im'][maxidx_l_cc]], label = r'Re$(\lambda)_{\rm max}$', s = 30, c = '#d95f02', marker = 'x', linewidth = 1.)
    except:
        maxidx_l = fig10_to_plot[s]['maxidx_l']
        ax2.scatter(fig10_to_plot[s]['lambdas_u_re'][maxidx_l], fig10_to_plot[s]['lambdas_u_im'][maxidx_l], label = r'Re$(\lambda)_{\rm max}$', s = 30, c = '#d95f02', marker = 'x', linewidth = 1.)
    ax2.set_xlabel(r'Re$(\lambda)$')
    ax2.set_ylabel(r'Im$(\lambda)$')

    if scheme == 'sbdf1':
        ax2.set_xlim(-3e-2, 1.2e-1)
        yb, yt = (-5e0, 5e0)
        ax2.set_ylim(1.1*yb, 1.1*yt)
        ax2.legend(loc = "lower right", borderaxespad=0.05, fontsize = 6, frameon=True)
    if scheme == 'sbdf2':
        ax2.set_xlim(-5e-6, 2e-5)
        ax2.set_xticks((0., 1e-5, 2e-5))
        ax2.set_xticklabels(('0','1e-5','2e-5'))
        yb, yt = (-1e1, 1e1)
        ax2.set_ylim(1.1*yb, 1.1*yt)
    if scheme == 'sbdf3':
        ax2.set_xlim(-6e0, 1.5e1)
        yb, yt = (-7e2, 7e2)
        ax2.set_ylim(1.1*yb, 1.1*yt)
    if scheme == 'sbdf4':
        ax2.set_xlim(-6e0, 6.5e1)
        yb, yt = (-9e2, 9e2)
        ax2.set_ylim(1.1*yb, 1.1*yt)
    if scheme == 'rk222':
        ax2.set_xlim(-2e-6, 1.3e-5)
        ax2.set_xticks((0., 5e-6, 1e-5))
        ax2.set_xticklabels(('0','5e-6','1e-5'))
        yb, yt = (-1e1, 1e1)
        ax2.set_ylim(1.1*yb, 1.1*yt)
    if scheme == 'rk443':
        ax2.set_xlim(-4e-6, 4e-6)
        ax2.set_xticks((-3e-6, 0., 3e-6))
        ax2.set_xticklabels(('-3e-6','0','3e-6'))
        yb, yt = (-3e1, 3e1)
        ax2.set_ylim(1.1*yb, 1.1*yt)
    
    ax2.axvline(x = 0., ymin = yb, ymax = yt, color = "black", linestyle = '--', zorder = -1)

plt.savefig('figure-10-vn-eigenspectra.eps')
plt.clf()
