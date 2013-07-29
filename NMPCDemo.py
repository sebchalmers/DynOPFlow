# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Runs an NMPC scheme on a simple power grid

"""

import DynOPFlow 
reload(DynOPFlow)
from DynOPFlow import *

################    DEFINE GRID TOPOLOGY    ##################
# Undirected connectivity graph: node i, node j, Zij
#Paper graph
NBus = 6 # Number of bus

#Connection pair (i,j), Impedance (undirected)
Graph = [  [ 0,1, 1+10j   ],
           [ 1,2, 1+10j   ],
           [ 1,3, 10+100j ],
           [ 3,5, 1+10j   ],
           [ 3,4, 1+10j   ]  ] 

#Define Net properties
Net = PowerGrid(NBus,Graph)
Net.Flow()

Net.PowerFlowBounds = {
                       'Vmin' :           [450.0 for k in range(NBus)],
                       'Vmax' :           [500.0 for k in range(NBus)],
                       'LineCurrentMax' : [50.0 for k in range(5)    ]
                      }


dt = 1.


#####  Define Hydro Plant #####
Hydro = Plant(States = ['h'], Inputs = ['qflow'], R = 0.1, Directionality = 'Bi', Bus = 4, label = 'Hydro')

etaT       =  0.8#1.25
etaP       =  0.75
A          =  1e-3
rho_air    =  1.2
rho_water  =  1e3
gravity    =  9.81
qTurbmax   =  2*6e-4   
qflow      =  Hydro.Inputs['qflow']
PP         =  Hydro.Inputs['Pcharge']
PT         =  Hydro.Inputs['Pdischarge']
PP_prev    =  Hydro.InputsPrev['Pcharge']
PT_prev    =  Hydro.InputsPrev['Pdischarge']
h          =  Hydro.States['h']

dh = (etaP*PP - PT/etaT)/(rho_water*gravity*A*h) + qflow/A 
Const = [PT/etaT - qTurbmax*rho_water*gravity*h]
Cost = (1/etaT - 1)*PT + (1 - etaP)*PP #+ 1e-1*(PP - PP_prev)**2  + 1e-1*(PT - PT_prev)**2

Hydro.setDynamics     (  RHS = dh, dt = dt    )
Hydro.setConstraints  (  Const                )
Hydro.setCost         (  Cost                 )

Net.addPlant(Hydro)

Hydro.LB['States','h'] = 5.
Hydro.UB['States','h'] = 20.
Hydro.UB['Inputs','Pcharge']     = 500.
Hydro.UB['Inputs','Pdischarge']  = 1000.



#####   Define Storage   #####  
Storage = Plant(States = ['E'], R = 0.1,  Directionality = 'Bi', Bus = 1, label = 'Storage')
etaC       = 0.9
etaD       = 0.95
tau        = 1e-6
Pcharge    = Storage.Inputs['Pcharge']
Pdischarge = Storage.Inputs['Pdischarge']
E          = Storage.States['E']

dEnergy = etaC*Pcharge - Pdischarge/etaD - tau*E
Storage.setDynamics( RHS = dEnergy, dt = dt )

Cost = (1/etaD - 1)*Pdischarge + (1 - etaC)*Pcharge #+ Pcharge*Pdischarge
Storage.setCost(Cost)
Net.addPlant(Storage)

Storage.LB['States','E']      = 0.
Storage.UB['States','E']      = 2e3
Storage.UB['Inputs','Pcharge']     = 250.
Storage.UB['Inputs','Pdischarge']  = 500.



#####  Define wind farm  #####  
Prated        = 1100 #Total rated power 
Wrated        = 10.
rho_air       = 1.23
A             = 2*Prated/(rho_air*0.47*Wrated**3)
CPmax         = .47
WindCurt      = 22
Tau           = 0.0
WindSpeedMean = 10.

Wind          = Plant(States = ['W'], Inputs = ['dW'], R = 0.1, Bus = 5, label = 'Wind')
PWind         = 0.5*rho_air*A*CPmax*Wind.States['W']**3

Const = []
Const.append(Wind.Inputs['Power'] - PWind)
Const.append(Wind.Inputs['Power']*(Wind.States['W']-WindCurt)/WindCurt/Prated - 1e-3)
Wind.setConstraints(Const)

#Wind random walk
dotW         = Wind.Inputs['dW'] - Tau*(Wind.States['W'] - WindSpeedMean)
Wind.setDynamics( RHS = dotW, dt = dt)
  
Net.addPlant(Wind)
Wind.UB['Inputs','Power']      = Prated


#####   Thermal   #####
ThermalRamp       = 200. 
Thermal           = Plant(Bus = 2, R = 0.1, label = 'Thermal')
ThermalPower      = Thermal.Inputs['Power']
ThermalPower_prev = Thermal.InputsPrev['Power']

Cost = 1e3*ThermalPower + (ThermalPower - ThermalPower_prev)**2
Thermal.setCost(Cost)

Const =   [    ThermalPower - ThermalPower_prev - ThermalRamp  ]  # ThermalPower - ThermalPower_prev <= ThermalRamp
Const.append( -ThermalPower + ThermalPower_prev - ThermalRamp  )  # - ThermalRamp <= ThermalPower - ThermalPower_prev
Thermal.setConstraints(Const)

Net.addPlant(Thermal)

Thermal.UB['Inputs','Power'] = 1000

#####   Load      ######
Load = Plant(Load = True, Bus = 0, label = 'Load')
Net.addPlant(Load)

Load.LB['Inputs',  'ActivePower'] = -1000
Load.LB['Inputs','ReactivePower'] =  -750
Load.UB['Inputs',  'ActivePower'] = -1000
Load.UB['Inputs','ReactivePower'] =  -750

# Impose current bounds on all plants
for plant in Net.PlantList:
    plant.UB['Inputs','CurrentReal'] =  5.
    plant.LB['Inputs','CurrentReal'] = -5.
    plant.UB['Inputs','CurrentImag'] =  5.
    plant.LB['Inputs','CurrentImag'] = -5.

#################    END OF NETWORK DEFINITION    ###########################

Horizon    = 48
Nsimulation = int(10*24)

Net.Profiles(Horizon + Nsimulation)

Nprofile = Net.Nprofile

dW = [rand.normalvariate(0,0.2) for k in range(Nprofile)]

LoadActivePower   = [300*np.cos(2*np.pi*k*dt/24.) - 1000  for k in range(Nprofile)]
LoadReactivePower = [0.75*LoadActivePower[k] for k in range(Nprofile)]

Net.Dispatch(Horizon = Horizon, Simulation = Nsimulation)

#Initial conditions (set inf in x0 to free the initial conditions)
u0 = Net.u0()
x0 = Net.x0()

u0['Thermal','Power']  = 0.
x0['Wind',   'W']      = 9.25
x0['Storage','E']      = 0.9*2e3
x0['Hydro',  'h']      = 0.9*20

#Make initial guess
init = Net.init()

init['States',:,'Wind','W'] = x0['Wind',   'W'] 
 
Net.LBInputProfiles['Hydro',:,'qflow'] = 6e-4
Net.UBInputProfiles['Hydro',:,'qflow'] = 6e-4


Net.LBInputProfiles['Wind',:,'dW']   = dW                                                                 
Net.UBInputProfiles['Wind',:,'dW']   = dW   
               
Net.LBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.LBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
Net.UBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.UBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
                                              
#Sol,_ = Net.DYNSolve(x0 = x0, u0 = u0, init = init)
#
#Net.ExtractInfo(Sol, PlantPower = True, BusPower = True, TotalPower = True)
#Net.DYNSolvePlot(Sol, dt = 1)
#
#assert(0==1)                                             
Traj, NMPC_Info = Net.NMPCSimulation(x0 = x0, u0 = u0, init = init, Simulation = Nsimulation) 
                       
#Plotting
Net.ExtractInfo(Traj)

Path = '/Users/sebastien/Desktop/Research/PowerFlowCodes/Paper/Figures/Simulations/Sim5'
SavedFigs = Net.DYNSolvePlot(Traj, dt = 1/24., Path = Path, LW = 2)          


## Create additional figures for the paper

    
    
 
