import numpy as np
import dedalus.public as d3
import h5py

scheme = 'sbdf2'
L = 10
alphas = np.logspace(np.log10(0.003), np.log10(0.02), 10) 
dts = np.logspace(np.log10(0.0006), np.log10(0.05), 22) 
i, j = (4, 10)
alpha = alphas[i]

bases = ['f', 'c']
Ns = np.array((128, 512))
colors = {}
for basis in bases:
    colors[basis] = []

data_dict = {}
for basis in bases:
    data_dict[basis] = {}
    for N in Ns:
        data_dict[basis][N] = {}
        if basis == 'f':
            output_name_f = f'{scheme}_output_f_{i}_{j}_{N}'
            file_in = h5py.File(output_name_f + '/' + output_name_f + '_s1.h5', 'r+')
            data_dict[basis][N]["ts_f"] = np.array(file_in['scales/sim_time'])
            data_dict[basis][N]["u_f"] = np.array(file_in['tasks']['u_f'])
            data_dict[basis][N]["max_u_f"] = np.max(np.array(file_in['tasks']['u_f']), axis = 1)
            data_dict[basis][N]["int_u_sq_f"] = np.array(file_in['tasks']['int_u_f_sq'])
        if basis == 'c':
            output_name_c = f'{scheme}_output_c_{i}_{j}_{N}'
            file_in = h5py.File(output_name_c + '/' + output_name_c + '_s1.h5', 'r+')
            data_dict[basis][N]["ts_c"] = np.array(file_in['scales/sim_time'])
            data_dict[basis][N]["u_c"] = np.array(file_in['tasks']['u_c'])
            data_dict[basis][N]["max_u_c"] = np.max(np.array(file_in['tasks']['u_c']), axis = 1)
            data_dict[basis][N]["int_u_sq_c"] = np.array(file_in['tasks']['int_u_c_sq'])


##### processing for ax1 #####
basis = 'f'
N = 512
ts = data_dict[basis][N]["ts_f"]
us = data_dict[basis][N]["u_f"]
cs = (1/3)*data_dict[basis][N]["max_u_f"]
# estimate number of domain crossings
deltats = np.diff(ts)
distance_estimate = np.sum(deltats * cs[:-1])
domain_crossings_estimate = distance_estimate / L
# find what the solution looked like after N crossings
idx0 = 0
ts_0 = ts[idx0]
us_0 = us[idx0,:]
idx300 = np.where(np.cumsum(deltats * cs[:-1])/L >= 300)[0][0] 
ts_300 = ts[idx300]
us_300 = us[idx300,:]
idx600 = np.where(np.cumsum(deltats * cs[:-1])/L >= 600)[0][0] 
ts_600 = ts[idx600]
us_600 = us[idx600,:]

#print("T_0", ts_0, "T_300", ts_300,"T_600", ts_600, "T_B", ts[-1])

# represent in complex fourier
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype = np.complex128)
xbasis = d3.ComplexFourier(xcoord, size = N, bounds = (-L/2, L/2))
x = dist.local_grid(xbasis)
k = xbasis.wavenumbers
u0 = dist.Field(name = 'u0', bases = xbasis)
u0['g'] = us_0
u300 = dist.Field(name = 'u300', bases = xbasis)
u300['g'] = us_300
u600 = dist.Field(name = 'u600', bases = xbasis)
u600['g'] = us_600
# shift to desired part of interval for visual
targets = [-L/3, 0, L/3]
x0 = x[np.argmax(u0['g'])]
p0 = targets[0] - x0
u0['c'] *= np.exp(-1j*k*p0)
x300 = x[np.argmax(u300['g'])]
p300 = targets[1] - x300
u300['c'] *= np.exp(-1j*k*p300)
x600 = x[np.argmax(u600['g'])]
p600 = targets[2] - x600
u600['c'] *= np.exp(-1j*k*p600)

# output processed data
ax1_to_plot = {}
ax1_to_plot['x'] = x
ax1_to_plot['u_T_0'] = u0['g']
ax1_to_plot['u_T_300'] = u300['g']
ax1_to_plot['u_T_600'] = u600['g']
np.save('fig1-ax1-processed.npy', ax1_to_plot)

##### processing for ax2 #####
ax2_to_plot = {}
for basis in bases:
    ax2_to_plot[basis] = {}
    for N in Ns: 
        ax2_to_plot[basis][N] = {}
        if basis == 'f':
            ax2_to_plot[basis][N]['t'] = data_dict[basis][N]["ts_f"]
            ax2_to_plot[basis][N]['L2'] = np.sqrt(data_dict[basis][N]["int_u_sq_f"])
        if basis == 'c':
            ax2_to_plot[basis][N]['t'] = data_dict[basis][N]["ts_c"]
            ax2_to_plot[basis][N]['L2'] = np.sqrt(data_dict[basis][N]["int_u_sq_c"])

# some values from ax1
ax2_to_plot['T_0'] = ts_0
ax2_to_plot['T_300'] = ts_300
ax2_to_plot['T_600'] = ts_600
basis = 'f'
N = 512
ax2_to_plot['L2_T_0'] = np.sqrt(data_dict[basis][N]["int_u_sq_f"][idx0])
ax2_to_plot['L2_T_300'] = np.sqrt(data_dict[basis][N]["int_u_sq_f"][idx300])
ax2_to_plot['L2_T_600'] = np.sqrt(data_dict[basis][N]["int_u_sq_f"][idx600])

np.save('fig1-ax2-processed.npy', ax2_to_plot)
