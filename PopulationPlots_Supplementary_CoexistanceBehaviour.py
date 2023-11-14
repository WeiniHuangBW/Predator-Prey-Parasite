# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:05:46 2022

@author: amymm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import odeint

# set rates
gx = 2
fy = 0.01
Scap = 0.0005
dx = 0.1
dy = 1
dz = 0.09
rx = 1
#rp = 1
re = 1
nz = 6
#Qx = 1
#Qy = 1
ky = 0.2
Kcap = 2000

# set Qx, Qy and rp rates
Qx = 0.8
Qy = 1
rp = 0.8

# set max time on plot
xmax = 100

# initial condition
a0 = [0,800,0,100,1000]

# data requirements
intervals = 100
domainsize = intervals * intervals * intervals
plotsize = intervals * intervals
order = 2
timeend = 500
threshold = 1

# time points
t = np.linspace(0,timeend, 10000)

# assign colours
IpreyColour = '#808080'
UpreyColour = '#00FFFF'
IpredatorColour = 'black'
UpredatorColour = '#FFFF00'
parasiteColour = '#FF00FF'

#set legend
legendline = [Line2D([0],[0], linewidth = 8,color = UpredatorColour),Line2D([0], [0], linewidth = 8, color = UpreyColour),Line2D([0],[0], linewidth = 8,color = IpredatorColour),Line2D([0], [0], linewidth = 8, color = IpreyColour),Line2D([0],[0], linewidth = 8,color = parasiteColour)]
legendnames = ['Uninfected predator', 'Uninfected prey', 'Infected predator', 'Infected prey', 'Parasite']

# data storage
outputFolder = r'C:/MyFiles/Masters/Masters/Project/POSTWORK/Code/DataCeci/Images'

outputPlot = outputFolder + '/populations_Qx_' + str(Qx) + '_Qy_' + str(Qy) + '_rp_' + str(rp) + '.png'
outputSimplex = outputFolder + '/simplex_Qx_' + str(Qx) + '_Qy_' + str(Qy) + '_rp_' + str(rp) + '.png'


#%%

# predator prey parasite ODE solver
def ODEsolve(a, t):
    # assign initial conditions
    xi = a[0]
    xu = a[1]
    yi = a[2]
    yu = a[3]
    z = a[4]
    # system of ODEs
    dxidt = -(xi + xu)*xi/Kcap - dx*xi + Qx*Scap*xu*z - fy* xi*(yi+yu)
    dxudt = gx*(rx*xi+xu) - (xi+xu)*xu/Kcap - dx*xu - Qx*Scap*xu*z - fy*xu*(yi+yu)
    dyidt = -dy*yi + Qy*fy*xi*yu
    dyudt = ky*fy*xu*(rp*yi+yu) + (re*ky*(1 - Qy) - (1 - rp*ky)*Qy)*fy*xi*yu + (rp*re*ky*(1 - Qy) + rp*rp*ky*Qy)*fy*xi*yi - dy*yu
    dzdt = -Qx*Scap*xu*z + nz*Qy*fy*xi*(yi + yu) - dz*z
    # output 
    dadt = [dxidt, dxudt, dyidt, dyudt, dzdt]
    return dadt




# solve ODE
z = odeint(ODEsolve, a0, t)

# create DataFrame of results
columnNames = ['xi', 'xu', 'yi', 'yu', 'z']
DFOutput = pd.DataFrame(z, columns = columnNames)
DFOutput['t'] = t

#%% Populations plot with parasite on alternate axis

# create figure    
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

# create alternate axis for parasite
ax2 = ax.twinx()

# plot populations over time
DFOutput.plot(x = 't', y = ['yu', 'xu','yi', 'xi'], ax = ax, color=[IpredatorColour,IpreyColour,UpredatorColour,UpreyColour], linewidth = 2, legend=False)#, columns={'xi_0', 'xu_0', 'yi_0', 'yu_0', 'z_0'})
DFOutput.plot(x = 't', y = ['z'], ax=ax2,color=[parasiteColour], linewidth = 2, legend=False)#, columns={'xi_0', 'xu_0', 'yi_0', 'yu_0', 'z_0'})

# plot kwargs
ax.set(xlim = (0, xmax), ylim = (0, 1.1*DFOutput[['xi', 'xu', 'yi', 'yu']].max(axis = 1).max(axis=0))) 
ax2.set(xlim = (0, xmax), ylim = (0, 1.1*DFOutput[['z']].max(axis = 1).max(axis=0))) 
ax.legend(legendline, legendnames, loc='upper center', prop={'size': 9.5}, ncol = 5, handleheight = 1, handlelength = 1.5)#, facecolor = "grey", labelcolor = "white")
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Host population size', fontsize=12)
ax2.set_ylabel('Parasite population size', color=parasiteColour, fontsize=12)
ax2.tick_params(axis='y', labelcolor=parasiteColour, colors = parasiteColour)
ax2.spines['right'].set_color(parasiteColour)
plt.savefig(outputPlot,dpi=300,bbox_inches='tight')
plt.show()
plt.rcParams.update(plt.rcParamsDefault)

#%% Simplex with scaling on parasite

# create figure
fig_sim = plt.figure()

# make copy of data
DF_popln_raw = pd.DataFrame(z)
DF_popln = DF_popln_raw.copy()

# evaluate total predator, prey and parasite populations
DF_popln['x'] = DF_popln[0] + DF_popln[1]
DF_popln['y'] = DF_popln[2] + DF_popln[3]
DF_popln['z'] = DF_popln[4]

# plot ternary
ax_sim = fig_sim.add_subplot(1, 1, 1, projection='ternary')
ax_sim.plot(DF_popln['x'], DF_popln['y'], DF_popln['z'], color='black')

# plot kwargs
ax_sim.grid()
ax_sim.set_tlabel('Prey')
ax_sim.set_llabel('Predator')
ax_sim.set_rlabel('Parasite', color = parasiteColour)
ax_sim.taxis.set_label_position('tick1')
ax_sim.laxis.set_label_position('tick1')
ax_sim.raxis.set_label_position('tick1')
ax_sim.taxis.set_tick_params(grid_color= UpreyColour)#colors=IpreyColour
ax_sim.laxis.set_tick_params(tick2On=True, grid_color= UpredatorColour)#, colors=IpredatorColour)
ax_sim.raxis.set_tick_params(tick2On=True, grid_color=parasiteColour, colors=parasiteColour)
ax_sim.spines['tside'].set_color(UpreyColour)
ax_sim.spines['lside'].set_color(UpredatorColour)
ax_sim.spines['rside'].set_color(parasiteColour)
plt.savefig(outputSimplex,dpi=300,bbox_inches='tight')
plt.show()



