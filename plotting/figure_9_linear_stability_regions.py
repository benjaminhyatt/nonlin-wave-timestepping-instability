"""
Plotting script for Figure 9 

Creates visuals of the linear stability regions for the IMEX
schemes considered in this work
"""
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

schemes = ['sbdf1', 'sbdf2', 'sbdf3', 'sbdf4', 'rk222', 'rk443']

def fmt(dt):
    return "{:.1}".format(dt)

def SBDF1(A1, A2, B):
    return np.abs((1 + 1j*B) / (1 - A1 - 1j*A2))

def SBDF2(A1, A2, B):
    A = A1 + 1j * A2
    p0 = 3/2 - A # highest order coeff
    p1 = -2 - 2*1j*B
    p2 = 1/2 + 1j*B
    Z = np.zeros((A2.shape[0], A2.shape[1])) # same as B.shape
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            p = np.array((p0[i, j], p1[i, j], p2[i, j]))
            Z[i, j] = np.max(np.abs(np.roots(p)))
    return Z

def SBDF3(A1, A2, B):
    A = A1 + 1j * A2
    p0 = 11/6 - A # highest order coeff
    p1 = -3 - 3*1j*B
    p2 = 3/2 + 3*1j*B
    p3 = -1/3 - 1j*B
    Z = np.zeros((A2.shape[0], A2.shape[0]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            p = np.array((p0[i, j], p1[i, j], p2[i, j], p3[i, j]))
            Z[i, j] = np.max(np.abs(np.roots(p)))
    return Z

def SBDF4(A1, A2, B):
    A = A1 + 1j * A2
    p0 = 25/12 - A # highest order coeff
    p1 = -4 - 4*1j*B
    p2 = 3 + 6*1j*B
    p3 = -4/3 - 4*1j*B
    p4 = 1/4 + 1j*B
    Z = np.zeros((A2.shape[0], A2.shape[0]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            p = np.array((p0[i, j], p1[i, j], p2[i, j], p3[i, j], p4[i, j]))
            Z[i, j] = np.max(np.abs(np.roots(p)))
    return Z

def RK222(A1, A2, B):
    A = A1 + 1j * A2
    gamma = (2 - np.sqrt(2))/2
    delta = 1 - 1/gamma/2
    KoverX0 = (1 + gamma*1j*B)*(1 - gamma*A)**(-1)
    X1overX0 = (1 + delta*1j*B + KoverX0*((1 - gamma)*A + (1 - delta)*1j*B))*(1 - gamma*A)**(-1)
    return np.abs(X1overX0)

def RK443(A1, A2, B):
    A = A1 + 1j * A2
    K1overX0 = (1 + 1j*B/2)*(1 - A/2)**(-1)
    K2overX0 = (1 + 1j*B*(11/18 + K1overX0/18) + (A/6)*K1overX0)*(1 - A/2)**(-1)
    K3overX0 = (1 + 1j*B*(5/6 - 5*K1overX0/6 + K2overX0/2) + A*(-K1overX0/2 + K2overX0/2))*(1 - A/2)**(-1)
    X1overX0 = (1 + 1j*B*(1/4 + 7*K1overX0/4 + 3*K2overX0/4 - 7*K3overX0/4) + A*(3*K1overX0/2 - 3*K2overX0/2 + K3overX0/2))*(1 - A/2)**(-1)
    return np.abs(X1overX0)

A1 = 0.
ctr_lvls = [0.5, 0.75, 1.]

##### make 6-panel figure #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 
t_mar, b_mar, l_mar, r_mar = (0.13, 0.26, 0.36, 0.13)

h_plot, w_plot = (1., 1.)
w_pad = 0.05 * w_plot
h_pad = w_pad

h_total = t_mar + h_plot + b_mar + h_plot + b_mar + h_plot + b_mar
w_total = l_mar + w_plot + w_pad + l_mar + w_plot + r_mar

fig_width = 5.5
scale = fig_width/w_total

fig = plt.figure(figsize = (scale * w_total, scale * h_total))

##### construct axs #####

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax1 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax2 = fig.add_axes([left, bottom, width, height])

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax3 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax4 = fig.add_axes([left, bottom, width, height])

left = (l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax5 = fig.add_axes([left, bottom, width, height])

left = (l_mar + w_plot + w_pad + l_mar) / w_total
bottom = 1 - (t_mar + h_plot + b_mar + h_plot + b_mar + h_plot) / h_total
width = w_plot / w_total
height = h_plot / h_total
ax6 = fig.add_axes([left, bottom, width, height])

axs = [ax1, ax2, ax3, ax4, ax5, ax6]

# plot axes
linestyles = ["dotted", "dotted", "solid"]
spacings = [(1, 2), (1, 1), (1, 0)] 

# ax1 - sbdf1
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
[X, Y] = np.meshgrid(x, y)

ax1.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax1.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = SBDF1(A1, X, Y)
    CS1 = ax1.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS1.collections:
        line.set_linestyle((0,spacings[i]))
ax1.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax1.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax1.set_xlim(x.min(), x.max())
ax1.set_ylim(y.min(), y.max())

labels = [r'$|\sigma| = 0.50$', r'$|\sigma| = 0.75$', r'$|\sigma| = 1$']
s0 = matplotlib.lines.Line2D([], [], color = "black", linestyle = linestyles[0], dashes = spacings[0], label = labels[0])
s1 = matplotlib.lines.Line2D([], [], color = "black", linestyle = linestyles[1], dashes = spacings[1], label = labels[1])
s2 = matplotlib.lines.Line2D([], [], color = "black", linestyle = linestyles[2], dashes = spacings[2], label = labels[2])
ax1.legend(handles = [s0, s1, s2], loc = "lower center", borderaxespad = 0.05, fontsize = 7, frameon = False)

# ax2 - sbdf2
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
[X, Y] = np.meshgrid(x, y)

ax2.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax2.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = SBDF2(A1, X, Y)
    CS2 = ax2.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS2.collections:
        line.set_linestyle((0,spacings[i]))
ax2.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax2.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax2.set_xlim(x.min(), x.max())
ax2.set_ylim(y.min(), y.max())

# ax3 - sbdf3
x = np.linspace(-20, 20, 200)
y = np.linspace(-20, 20, 200)
[X, Y] = np.meshgrid(x, y)

ax3.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax3.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = SBDF3(A1, X, Y)
    CS3 = ax3.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS3.collections:
        line.set_linestyle((0,spacings[i]))
ax3.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax3.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax3.set_xlim(x.min(), x.max())
ax3.set_ylim(y.min(), y.max())

# ax4 - sbdf4
x = np.linspace(-20, 20, 200)
y = np.linspace(-20, 20, 200)
[X, Y] = np.meshgrid(x, y)
ax4.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax4.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = SBDF4(A1, X, Y)
    CS4 = ax4.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS4.collections:
       line.set_linestyle((0,spacings[i]))
ax4.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax4.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax4.set_xlim(x.min(), x.max())
ax4.set_ylim(y.min(), y.max())

# ax5 - rk222
x = np.linspace(-20, 20, 200)
y = np.linspace(-20, 20, 200)
[X, Y] = np.meshgrid(x, y)
ax5.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax5.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = RK222(A1, X, Y)
    CS5 = ax5.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS5.collections:
        line.set_linestyle((0,spacings[i]))
ax5.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax5.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax5.set_xlim(x.min(), x.max())
ax5.set_ylim(y.min(), y.max())

# ax6 - rk443
x = np.linspace(-20, 20, 200)
y = np.linspace(-20, 20, 200)
[X, Y] = np.meshgrid(x, y)

ax6.hlines(0, x.min(), x.max(), linestyle = "dashed", color = "grey", lw = 0.4)
ax6.vlines(0, y.min(), y.max(), linestyle = "dashed", color = "grey", lw = 0.4)
for i, lvl in enumerate(ctr_lvls):
    Z = RK443(A1, X, Y)
    CS6 = ax6.contour(X, Y, Z, levels = np.array([lvl]), colors = "black", linewidths = 0.7, linestyles = linestyles[i])
    for line in CS6.collections:
        line.set_linestyle((0,spacings[i]))
ax6.set_xlabel(r"${\rm Im}(a_{\rm IM})$")
ax6.set_ylabel(r"${\rm Im}(a_{\rm EX})$")
ax6.set_xlim(x.min(), x.max())
ax6.set_ylim(y.min(), y.max())

for sch, ax in enumerate(axs):
    scheme = schemes[sch]
    ha = 'right'
    xx = 0.9
    ax.text(xx, 0.95, schemes[sch].upper(), ha = ha, va = 'top', transform=ax.transAxes)

plt.savefig(f"figure-9-stability-regions.eps")
plt.clf()

