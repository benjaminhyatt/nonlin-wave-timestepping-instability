import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import cmasher as cmr 

L = 10
bases = ['f', 'c']
colors = {}
for basis in bases:
    colors[basis] = []

basis = 'f'
colors[basis].append(matplotlib.colors.to_rgba('#80cdc1'))
colors[basis].append(matplotlib.colors.to_rgba('#018571'))
basis = 'c'
colors[basis].append(matplotlib.colors.to_rgba('#dfc27d'))
colors[basis].append(matplotlib.colors.to_rgba('#a6611a'))

##### load in processed data for plotting ##### 
ax1_to_plot = np.load('fig1-ax1-processed.npy', allow_pickle = True)[()]
ax2_to_plot = np.load('fig1-ax2-processed.npy', allow_pickle = True)[()]

##### construct fig #####
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.size'] = 8 
plt.rcParams['figure.dpi'] = 600 
fontsize = 8 

t_mar, b_mar, l_mar, r_mar = (0.30, 0.36, 0.44, 0.12)
golden_mean = (np.sqrt(5) - 1.) / 2.
h_plot, w_plot = (1., 1. / golden_mean)
w_pad = l_mar

h_total = t_mar + h_plot + b_mar
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

axs = [ax1, ax2]

##### ax1 #####
ax1.spines[['right', 'top']].set_visible(False)
ax1.plot(ax1_to_plot['x'], ax1_to_plot['u_T_0'], color = matplotlib.colors.to_rgba('#6baed6'), lw = 1, label = r"$t = T_0$")
ax1.plot(ax1_to_plot['x'], ax1_to_plot['u_T_300'], color = matplotlib.colors.to_rgba('#3182bd'), lw = 1, label = r"$t = T_{300}$")
ax1.plot(ax1_to_plot['x'], ax1_to_plot['u_T_600'], color = matplotlib.colors.to_rgba('#08519c'), lw = 1, label = r"$t = T_{600}$")

ax1.legend(loc = 'upper left', borderaxespad=0.05, fontsize = fontsize, frameon=False)

ax1.set_xlim(-L/2, L/2)
ax1.set_ylim(0., 3.)
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$u$")
ax1.set_xticks(())
ax1.set_yticks((0, 1.5, 3.))


##### ax2 #####
basis = 'f'
n, N = (0, 128)
ax2.plot(ax2_to_plot[basis][N]['t'], ax2_to_plot[basis][N]['L2'], color = colors[basis][n], lw = 1, label = rf"$N_F$ = {N}")

basis = 'f'
n, N = (1, 512)
ax2.plot(ax2_to_plot[basis][N]['t'], ax2_to_plot[basis][N]['L2'], color = colors[basis][n], lw = 1, linestyle = "dashed", dashes=(5, 5), label = rf"$N_F$ = {N}")

basis = 'c'
n, N = (1, 512)
ax2.plot(ax2_to_plot[basis][N]['t'], ax2_to_plot[basis][N]['L2'], color = colors[basis][n], lw = 1, linestyle = "dashed", dashes=(2.5, 5), label = rf"$N_C$ = {N}")

ax2.scatter(ax2_to_plot['T_0'], ax2_to_plot['L2_T_0'], s = 8, color = matplotlib.colors.to_rgba('#6baed6'), marker = "o")
ax2.annotate(r"$T_0$", (ax2_to_plot['T_0'], ax2_to_plot['L2_T_0']), textcoords="offset points", xytext=(0,4), ha='center')
ax2.vlines(ax2_to_plot['T_0'], 0, ax2_to_plot['L2_T_0'], color = matplotlib.colors.to_rgba('#6baed6'), lw = 0.5, ls = "dashed")
ax2.scatter(ax2_to_plot['T_300'], ax2_to_plot['L2_T_300'], s = 8, color = matplotlib.colors.to_rgba('#3182bd'), marker = "o")
ax2.annotate(r"$T_{300}$", (ax2_to_plot['T_300'], ax2_to_plot['L2_T_300']), textcoords="offset points", xytext=(0,5), ha='center')
ax2.vlines(ax2_to_plot['T_300'], 0, ax2_to_plot['L2_T_300'], color = matplotlib.colors.to_rgba('#3182bd'), lw = 0.5, ls = "dashed")
ax2.scatter(ax2_to_plot['T_600'], ax2_to_plot['L2_T_600'], s = 8, color = matplotlib.colors.to_rgba('#08519c'), marker = "o")
ax2.annotate(r"$T_{600}$", (ax2_to_plot['T_600'], ax2_to_plot['L2_T_600']), textcoords="offset points", xytext=(-10,4), ha='center')
ax2.vlines(ax2_to_plot['T_600'], 0, ax2_to_plot['L2_T_600'], color = matplotlib.colors.to_rgba('#08519c'), lw = 0.5, ls = "dashed")

ax2.set_ylim(5e-1, 5e1)
ax2.set_yscale("log")

ax2.xaxis.set_major_locator(matplotlib.ticker.FixedLocator((0, 3e3, 6e3, 9e3)))
ax2.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(n=3))

ax2.set_xlabel(r"$t$")
ax2.set_ylabel(r"$L_2(u)$")
ax2.legend(loc = 'upper left', borderaxespad=0.05, fontsize = fontsize, frameon=False)

plt.savefig(f"figure-1-instability-visual.eps")
plt.clf()



