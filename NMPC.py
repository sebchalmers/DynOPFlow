# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Runs a Monte Carlo simulation of NMPC loops with random wind profiles, various horizon length, and WITH forecast uncertainty


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
Nsimulation = int(10*24)

#####  Define Hydro Plant #####
Hydro = Plant(States = ['WaterHeight'], Inputs = ['qflow'], R = 0.1, Directionality = 'Bi', Bus = 4, label = 'Hydro')

etaT      =  1.25
etaP      =  0.75
A         =  1e-3
rho_air   =  1.2
rho_water =  1e3
gravity   =  9.81
qTurbmax  =  2*6e-4   
qflow     =  Hydro.Inputs['qflow']
PPump     =  Hydro.Inputs['Pcharge']
PTurb     =  Hydro.Inputs['Pdischarge']
h         =  Hydro.States['WaterHeight']

dh = (etaP*PPump - etaT*PTurb)/(rho_water*gravity*A*h) + qflow/A 
Const = [etaT*PTurb - qTurbmax*rho_water*gravity*h]
Cost = (etaT - 1)*PTurb + (1 - etaP)*PPump + PPump*PTurb

Hydro.setDynamics     (  RHS = dh, dt = dt    )
Hydro.setConstraints  (  Const                )
Hydro.setCost         (  Cost                 )

#TConst = (15-h)
#Hydro.setConstraints(TConst, Terminal = True)
#TCost = -100*h
#Hydro.setCost(TCost, Terminal = True)

Hydro.addPlant(Net)

Hydro.LB['States','WaterHeight'] = 5.
Hydro.UB['States','WaterHeight'] = 20.
Hydro.UB['Inputs','Pcharge']     = 500.
Hydro.UB['Inputs','Pdischarge']  = 1000.

#Hydro.LB['Inputs','qflow']  = 6e-4
#Hydro.UB['Inputs','qflow']  = 6e-4


#####   Define Storage   #####  
Storage = Plant(States = ['Energy'], R = 0.1,  Directionality = 'Bi', Bus = 1, label = 'Storage')
etaC       = 0.9
etaD       = 1.1
tau        = 1e-6
Pcharge    = Storage.Inputs['Pcharge']
Pdischarge = Storage.Inputs['Pdischarge']
E          = Storage.States['Energy']

dEnergy = etaC*Pcharge - etaD*Pdischarge - tau*E
Storage.setDynamics( RHS = dEnergy, dt = dt )


Cost = (etaD - 1)*Pdischarge + (1 - etaC)*Pcharge + Pcharge*Pdischarge
Storage.setCost(Cost)

Storage.addPlant(Net)

Storage.LB['States','Energy']      = 0.
Storage.UB['States','Energy']      = 2e3
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

Wind          = Plant(States = ['WindSpeed'], Inputs = ['dWindSpeed'], R = 0.1, Bus = 5, label = 'Wind')

#Wind availaibility
PWind         = 0.5*rho_air*A*CPmax*Wind.States['WindSpeed']**3

Const = []
Const.append(Wind.Inputs['Power'] - PWind)
Const.append(Wind.Inputs['Power']*(Wind.States['WindSpeed']-WindCurt)/WindCurt/Prated - 1e-3)
Wind.setConstraints(Const)

#Wind random walk
dWind         = Wind.Inputs['dWindSpeed'] - Tau*(Wind.States['WindSpeed'] - WindSpeedMean)
Wind.setDynamics( RHS = dWind, dt = dt)
            
Wind.addPlant(Net)

Wind.UB['Inputs','Power']      = Prated
#Wind.LB['Inputs','dWindSpeed'] = 0
#Wind.UB['Inputs','dWindSpeed'] = 0

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



Thermal.addPlant(Net)

Thermal.UB['Inputs','Power'] = 1000

#####   Load      ######
Load = Plant(Load = True, Bus = 0, label = 'Load')
Load.addPlant(Net)
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
Nsimulation = 10*24

Net.Profiles(Horizon + Nsimulation)

Nprofile = Net.Nprofile

dWind             = []
LoadActivePower   = []
LoadReactivePower = []


dWind = [rand.normalvariate(0,0.0) for k in range(Nprofile)]
LoadMean = [-300*np.cos(2*np.pi*k*dt/24.) + 1000  for k in range(Nprofile)]
#LoadMean = [0*np.cos(2*np.pi*k*dt/24.) + 1000  for k in range(Nprofile)]

LoadStd  =  0.0
dWindStd =  0.0

LoadActivePower   = [min(0,rand.normalvariate(-LoadMean[k],LoadStd)) for k in range(Nprofile)]
LoadReactivePower = [0.75*LoadActivePower[k] for k in range(Nprofile)]




Net.Dispatch(Horizon = Horizon, Simulation = Nsimulation)

#Initial conditions (set inf in x0 to free the initial conditions)
u0 = Net.u0()
x0 = Net.x0()

u0['Thermal','Power']         = 0.
x0['Wind',   'WindSpeed']     = 9.25
x0['Storage','Energy']        = 0.9*2e3
x0['Hydro',  'WaterHeight']   = 0.9*20

#Make initial guess
init = Net.init()

init['States',:,'Wind','WindSpeed'] = 10.  
 
Net.LBInputProfiles['Hydro',:,'qflow'] = 6e-4
Net.UBInputProfiles['Hydro',:,'qflow'] = 6e-4


Net.LBInputProfiles['Wind',:,'dWindSpeed']   = [dWind[k] + rand.normalvariate(0,dWindStd)   for k in range(Nprofile)]                                                                     
Net.UBInputProfiles['Wind',:,'dWindSpeed']   = [dWind[k] + rand.normalvariate(0,dWindStd)   for k in range(Nprofile)]      
               
Net.LBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.LBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
Net.UBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.UBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
                                              
#Solve initial guess

Sol, stats = Net.DYNSolve(x0 = x0, u0 = u0, time = 0, init = init)
init = Sol

Net.ExtractInfo(Sol, PlantPower = True, BusPower = True, TotalPower = True)
#Net.DYNSolvePlot(Sol, dt = 1)


Traj, NMPC_Info = Net.NMPCSimulation(x0 = x0, u0 = u0, init = init, Simulation = Nsimulation) 
                       

#Exctract info
Net.ExtractInfo(Traj, PlantPower = True, BusPower = True, TotalPower = True)
#Net.DYNSolvePlot(Traj, dt = 1)          




    
    
