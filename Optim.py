# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Runs an NMPC scheme on a simple power grid

"""
import os as os
import numpy as np
import sys
sys.path.append('/Users/sebastien/Desktop/DynOPFlow')
import DynOPFlow2
reload(DynOPFlow2)
from DynOPFlow2 import *

#from scipy import linalg
def null(A, eps=1e-15):
    u, s, vh = np.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)

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
                       'Vmin' :           [10.0 for k in range(NBus)],
                       'Vmax' :           [500.0 for k in range(NBus)],
                       'LineCurrentMax' : [5. for k in range(5)    ]
                      }


dt = 1.


#####  Define Hydro Plant #####
#Hydro = Plant(States = ['h'], Inputs = ['qflow'], R = 0.1, Directionality = 'Bi', Bus = 4, label = 'Hydro')
Hydro = Plant(States = ['h'],  R = 0.1, Directionality = 'Bi', Bus = 4, label = 'Hydro')

etaT       =  0.8
etaP       =  0.75
A          =  1e-3
rho_air    =  1.2
rho_water  =  1e3
gravity    =  9.81
qTurbmax   =  2*6e-4   
qflow      =  6e-4#Hydro.Inputs['qflow']
PP         =  Hydro.Inputs['Pcharge']
PT         =  Hydro.Inputs['Pdischarge']
PP_prev    =  Hydro.InputsPrev['Pcharge']
PT_prev    =  Hydro.InputsPrev['Pdischarge']
h          =  Hydro.States['h']

dh = (etaP*PP - PT/etaT)/(rho_water*gravity*A*h) + qflow/A 
Const = [PT/etaT - qTurbmax*rho_water*gravity*h]
Cost = (1/etaT - 1)*PT + (1 - etaP)*PP 

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

Cost = (1/etaD - 1)*Pdischarge + (1 - etaC)*Pcharge 
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
#for plant in Net.PlantList:
#    plant.UB['Inputs','CurrentReal'] =  5
#    plant.LB['Inputs','CurrentReal'] = -5
#    plant.UB['Inputs','CurrentImag'] =  5
#    plant.LB['Inputs','CurrentImag'] = -5

#################    END OF NETWORK DEFINITION    ###########################
def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)
        
Horizon    = 24
#Nsimulation = int(10*24)

Net.Profiles(Horizon)

Nprofile = Net.Nprofile

dW = [rand.normalvariate(0,0.0) for k in range(Nprofile)]

LoadActivePower   = [300*np.cos(2*np.pi*k*dt/24.) - 1000  for k in range(Nprofile)]
LoadReactivePower = [0.75*LoadActivePower[k] for k in range(Nprofile)]

Net.Dispatch(Horizon = Horizon),# Simulation = Nsimulation)

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
 
#Net.LBInputProfiles['Hydro',:,'qflow'] = 6e-4
#Net.UBInputProfiles['Hydro',:,'qflow'] = 6e-4

Net.LBInputProfiles['Wind',:,'dW']   = dW                                                                 
Net.UBInputProfiles['Wind',:,'dW']   = dW   
               
Net.LBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.LBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
Net.UBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.UBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
                                              
Sol,_ = Net.DYNSolve(x0 = x0, u0 = u0, init = init)

Net._HessOptDispatch.setInput(Net.OptDispatch.output('x'),0)
Net._HessOptDispatch.setInput(1.,1)
Net._HessOptDispatch.setInput(1.,2)
Net._HessOptDispatch.setInput(Net.OptDispatch.output('lam_g'),3)
Net._HessOptDispatch.evaluate()
H = Net._HessOptDispatch.output()

Net._JacOptDispatch.setInput(Net.OptDispatch.output('x'),0)
Net._JacOptDispatch.evaluate()
J = Net._JacOptDispatch.output()

#g_sol = Net.gOptDispatch(Net.OptDispatch.output('g'))

g_sol = Net.OptDispatch.output('g')
V_sol = Net.VOptDispatch.cat
lam_g_sol = Net.OptDispatch.output('lam_g')
i_active = []
g_active = []
for i in range(g_sol.shape[0]):
    if (Net.ubg.cat[i]-g_sol[i] < lam_g_sol[i]) or (g_sol[i]-Net.lbg.cat[i] < lam_g_sol[i]):
        i_active.append(i)
        if (Net.ubg.cat[i]-g_sol[i] < lam_g_sol[i]):
            g_active.append(Net.ubg.cat[i]-g_sol[i])
        else:
            g_active.append(g_sol[i]-Net.lbg.cat[i])   


i_bound = []
Jbound = []
for i in range(Sol.shape[0]):
    if (Sol.cat[i] - Net.lbV.cat[i]  < 1e-10) or (Net.ubV.cat[i] - Sol.cat[i] < 1e-10):
        i_bound.append(i)
        if (Sol.cat[i] - Net.lbV.cat[i]  < 1e-10):
            Newline = np.zeros([1,Sol.shape[0]])
            Newline[0,i] = 1.
            Jbound.append(Newline)
            g_active.append(Sol.cat[i] - Net.lbV.cat[i])
        else:
            Newline = np.zeros([1,Sol.shape[0]])
            Newline[0,i] = -1.
            Jbound.append(Newline)
            g_active.append(Net.ubV.cat[i] - Sol.cat[i])
        
Jbound = np.concatenate(Jbound,axis = 0)
g_active = np.concatenate(g_active,axis=0)
                         
iJ = i_active+i_bound

Net._HessOptDispatch.setInput(Net.OptDispatch.output('x'),0)
#Net._HessOptDispatch.setInput(Net.OptDispatch.output('x'),0)
#Net._HessOptDispatch.setInput(Net.OptDispatch.output('x'),0)
Net._HessOptDispatch.setInput(Net.OptDispatch.output('lam_g'),3)
Net._HessOptDispatch.evaluate()
H = Net._HessOptDispatch.output()

Net._JacOptDispatch.setInput(Net.OptDispatch.output('x'),0)
Net._JacOptDispatch.evaluate()
J = Net._JacOptDispatch.output()


J_active = [np.array(J[i_active,:])]
J_active.append(Jbound)
J_active = np.concatenate(J_active,axis=0)

KKT = np.zeros([J_active.shape[1] + J_active.shape[0],J_active.shape[1] + J_active.shape[0]])
KKT[:J_active.shape[1],:J_active.shape[1]] = H+np.eye(J_active.shape[1])
KKT[J_active.shape[1]:,:J_active.shape[1]] = J_active
KKT[:J_active.shape[1],J_active.shape[1]:] = J_active.T

Net._JacCostOptDispatch.setInput(Net.OptDispatch.output('x'),0)
Net._JacCostOptDispatch.evaluate()
JCost = Net._JacCostOptDispatch.output()

RHS = np.zeros([J_active.shape[1] + J_active.shape[0],1])
RHS[:J_active.shape[1],:] = JCost
RHS[J_active.shape[1]:,:] = g_active

dX = np.linalg.solve(KKT,RHS)[:J_active.shape[1]]
print "Norm 2 of dX:", np.sqrt(mul(dX.T,dX))


E,W = np.linalg.eig(KKT)
E = np.real(E)
print "Min eigenvalue of KKT mat", np.min(np.abs(E))

U,S,V = np.linalg.svd(J_active)
print "Min Singular value of J mat", np.min(S)


NJ = null(J_active.T,eps=1e-2)
#print "Null space: ", NJ

ising = []
for col in [0]:
    for k in range(len(NJ[:,col])):
        if (np.abs(NJ[k,col]) > 1e-1):
            ising.append(k)
        

Jsing = J_active[ising,:]
plt.subplot(2,1,1)
plt.spy(Jsing)
plt.subplot(2,1,2)
plt.spy(J)
plt.show()



    
    
 
