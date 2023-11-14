# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:05:46 2022

@author: amymm
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:05:46 2022

@author: amymm
"""

import numpy as np
import pandas as pd
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

# initial condition
a0 = [0,800,0,100,1000]

# data requirements
intervals = 100
plotsize = intervals * intervals
domainsize = intervals * intervals * intervals
order = 2
timeend = 100
threshold = 1

# time points
t = np.linspace(0,timeend, 1000)

# data storage
outputFolder = r'/output/'

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

#%%

# create data storage
dfOutput = pd.DataFrame()

# solve ODE for each element in the domain
for i in range (0, domainsize):

    # create parameter space values
    Qx = round((int(np.floor(i / (intervals * intervals))) + 1) * (1 / intervals), order)
    Qy = round((int(np.floor(i / intervals)) % intervals + 1) * (1 / intervals), order)
    rp = round((i % intervals + 1) * (1 / intervals), order)
    
    # solve ODE
    z = odeint(ODEsolve, a0, t)
    
    # Determine if its a steady state
    dfTemp = pd.DataFrame(z[-250:]).agg([np.mean, np.std, np.min, np.max]).transpose()
    dfTemp['CoV'] = dfTemp['std'] / dfTemp['mean']
    dfTemp['SteadyState']=np.where((dfTemp['CoV'] >= 0.1) & (dfTemp['amax']-dfTemp['amin'] > 1), 0, 1)
    steadyState = dfTemp['SteadyState'].min()
    
    # Determine if any populations are zero
    prey_i = round(z[-1,0])
    prey_u = round(z[-1,1])
    prey_t = prey_i + prey_u
    
    pred_i = round(z[-1,2])
    pred_u = round(z[-1,3])
    pred_t = pred_i + pred_u
    
    parasite = round(z[-1,4])
    
    ## CLASSIFICATIONS
    if steadyState == 0:
        # Divergent
        speciesExist = 4
    elif prey_t > 0 and pred_t > 0 and parasite > 0:
        # Pred + Prey + Parasite
        speciesExist = 1
    elif prey_t > 0 and pred_t > 0 and parasite == 0: 
        # Pred + Prey 
        speciesExist = 2
    elif prey_t == 0 and pred_t > 0 and parasite >= 0: 
        # Predator Only
        speciesExist = 3
    elif prey_t > 0 and pred_t == 0 and parasite >= 0: 
        # Prey only
        speciesExist = 99
    else:
        # Check for other cases
        speciesExist = 999
    
    ## SubType Dominance
    if speciesExist != 1:
        subtypeDominance = 5
    elif prey_u >= prey_i and pred_u >= pred_i:
        subtypeDominance = 1
    elif prey_u >= prey_i and pred_u < pred_i:
        subtypeDominance = 3
    elif prey_u < prey_i and pred_u >= pred_i:
        subtypeDominance = 2
    elif prey_u < prey_i and pred_u < pred_i:
        subtypeDominance = 4
    else:
        subtypeDominance = 99
    
    # Combine into a data frame
    results_df = pd.DataFrame({'Qx': [Qx]
                                 , 'Qy': [Qy]
                                 , 'rp': [rp]
                                 , 'xi': [prey_u]
                                 , 'xu': [prey_i]
                                 , 'yi': [pred_u]
                                 , 'yu': [pred_i]
                                 , 'x': [prey_t]
                                 , 'y': [pred_t]
                                 , 'z': [parasite]
                                 , 'speciesExist': [speciesExist]
                                 , 'subtypeDominance': [subtypeDominance]})
                                
    dfOutput = pd.concat([dfOutput, results_df])
    if (i % 1000 == 0):
        print('iteration: ', i)

# store data
dfOutput.to_csv(outputFolder + 'ParameterSpaceData.csv')

#%% create plot data for Qx_Qy, Qx_rp and Qy_rp

# copy parameter space data
dfOutput_reduced = dfOutput[['Qx', 'Qy', 'rp', 'xi', 'xu', 'yi', 'yu', 'z', 'speciesExist', 'subtypeDominance']].copy()

# create temp data storage
data_Qx_Qy_transpose = pd.DataFrame()
data_Qx_rp_transpose = pd.DataFrame()
data_Qy_rp_transpose = pd.DataFrame()

# group data by Qx_Qy, Qx_rp and Qy_rp
for i in range (0, domainsize):
    if (dfOutput['rp'][i] == 1):
        data_Qx_Qy_transpose = pd.concat([data_Qx_Qy_transpose, dfOutput_reduced.iloc[i]], axis = 1)
    if (dfOutput['Qy'][i] == 1):
        data_Qx_rp_transpose = pd.concat([data_Qx_rp_transpose, dfOutput_reduced.iloc[i]], axis = 1)
    if (dfOutput['Qx'][i] == 1):
        data_Qy_rp_transpose = pd.concat([data_Qy_rp_transpose, dfOutput_reduced.iloc[i]], axis = 1)
    if (i % 100 == 0):
        print('data iteration: ', i)

# create Qx_Qy, Qx_rp and Qy_rp data
data_Qx_Qy = data_Qx_Qy_transpose.T.copy()
data_Qx_rp = data_Qx_rp_transpose.T.copy()
data_Qy_rp = data_Qy_rp_transpose.T.copy()

# store data
data_Qx_Qy['Qx', 'Qy', 'speciesExist'].to_csv(outputFolder + 'Figure_3_Qx_Qy.png')
data_Qx_Qy['Qx', 'Qy', 'subtypeDominance'].to_csv(outputFolder + 'Figure_4_Qx_Qy.png')
data_Qx_rp['Qx', 'rp', 'speciesExist'].to_csv(outputFolder + 'Figure_3_Qx_rp.png')
data_Qx_rp['Qx', 'rp', 'subtypeDominance'].to_csv(outputFolder + 'Figure_4_Qx_rp.png')
data_Qy_rp['Qy', 'rp', 'speciesExist'].to_csv(outputFolder + 'Figure_3_Qy_rp.png')
data_Qy_rp['Qy', 'rp', 'subtypeDominance'].to_csv(outputFolder + 'Figure_4_Qy_rp.png')


