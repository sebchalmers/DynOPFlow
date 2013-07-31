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
#import DynOPFlow
#reload(DynOPFlow)
from DynOPFlow import *

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

#Graph = [  [ 0,1, 2+15j   ],
#           [ 1,2, 3+5j   ],
#           [ 1,3, 12+93j ],
#           [ 3,5, 1.5+13j   ],
#           [ 3,4, 0.5+12j   ]  ] 

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
#Hydro = Plant(States = ['h'], Inputs = ['qflow'], R = 0.0, Directionality = 'Bi', Bus = 4, label = 'Hydro')
Hydro = []
HydroLabels = []
for i in [4]:
                       HydroLabels.append('Hydro'+str(i))
                       Hydro.append(Plant(States = ['h'],  R = 0.0, Directionality = 'Bi', Bus = i, label = HydroLabels[-1]))                       
                       etaT       =  0.8
                       etaP       =  0.75
                       A          =  1e-3
                       rho_air    =  1.2
                       rho_water  =  1e3
                       gravity    =  9.81
                       qTurbmax   =  2*6e-4   
                       qflow      =  6e-4#Hydro.Inputs['qflow']
                       PP         =  Hydro[-1].Inputs['Pcharge']
                       PT         =  Hydro[-1].Inputs['Pdischarge']
                       PP_prev    =  Hydro[-1].InputsPrev['Pcharge']
                       PT_prev    =  Hydro[-1].InputsPrev['Pdischarge']
                       h          =  Hydro[-1].States['h']
                       
                       dh = (etaP*PP - PT/etaT)/(rho_water*gravity*A*h) + qflow/A 
                       Const = [PT/etaT - qTurbmax*rho_water*gravity*h]
                       Cost = (1/etaT - 1)*PT + (1 - etaP)*PP 
                       
                       Hydro[-1].setDynamics     (  RHS = dh, dt = dt    )
                       Hydro[-1].setConstraints  (  Const                )
                       Hydro[-1].setCost         (  Cost                 )
                       
                       Net.addPlant(Hydro[-1])
                       
                       Hydro[-1].LB['States','h'] = 5.
                       Hydro[-1].UB['States','h'] = 20.
                       Hydro[-1].UB['Inputs','Pcharge']     = 500.
                       Hydro[-1].UB['Inputs','Pdischarge']  = 1000.

#Net.addPlant(Hydro)

#####   Define Storage   #####  
Storage = Plant(States = ['E'], R = 0.0,  Directionality = 'Bi', Bus = 1, label = 'Storage')
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

Wind          = Plant(States = ['W'], Inputs = ['dW'], R = 0.0, Bus = 5, label = 'Wind')
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
Thermal = []
ThermalLabels = []
for k in [2]:
                       ThermalLabels.append('Thermal'+str(k))          
                       ThermalRamp       = 200. 
                       Thermal.append(      Plant(Bus = k, R = 0.0, label = ThermalLabels[-1]))
                       ThermalPower      = Thermal[-1].Inputs['Power']
                       ThermalPower_prev = Thermal[-1].InputsPrev['Power']
                       
                       Cost = (ThermalPower - ThermalPower_prev)**2 + 1e1*ThermalPower
                       Thermal[-1].setCost(Cost)
                       
                       Const =   [    ThermalPower - ThermalPower_prev - ThermalRamp  ]  # ThermalPower - ThermalPower_prev <= ThermalRamp
                       Const.append( -ThermalPower + ThermalPower_prev - ThermalRamp  )  # - ThermalRamp <= ThermalPower - ThermalPower_prev
                       Thermal[-1].setConstraints(Const)
                       
                       Net.addPlant(Thermal[-1])
                       
                       Thermal[-1].UB['Inputs','Power'] = 1000           

#####   Load      ######
Load = Plant(Load = True, Bus = 0, label = 'Load')
#ActivePower   = Load.Inputs[  'ActivePower']
#ReactivePower = Load.Inputs['ReactivePower']
#Cost = 1e6*(ActivePower + 1300)**2 + 1e6*(ActivePower + 950)**2
#Load.setCost(Cost)

Net.addPlant(Load)


Load.LB['Inputs',  'ActivePower'] = -1300
Load.LB['Inputs','ReactivePower'] =  -950
Load.UB['Inputs',  'ActivePower'] = -1300
Load.UB['Inputs','ReactivePower'] =  -950

#Load.LB['Inputs',  'ActivePower'] = -1000
#Load.LB['Inputs','ReactivePower'] =  -750
#Load.UB['Inputs',  'ActivePower'] = -1000
#Load.UB['Inputs','ReactivePower'] =  -750

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

for key in ThermalLabels:
    u0[key,'Power']  = 0.
x0['Wind',   'W']      = 9.25
x0['Storage','E']      = 0.9*2e3

for key in HydroLabels:
    x0[key,  'h']      = 0.9*20
#x0['Hydro2',  'h']      = 0.9*20

#Make initial guess
init = Net.init()

init['States',:,'Wind','W'] = x0['Wind',   'W'] 
 


Net.LBProfiles['Inputs',:,'Wind','dW']   = dW                                                                
Net.UBProfiles['Inputs',:,'Wind','dW']   = dW   
               
Net.LBProfiles['Inputs',:,'Load','ActivePower']   = LoadActivePower
Net.LBProfiles['Inputs',:,'Load','ReactivePower'] = LoadReactivePower
Net.UBProfiles['Inputs',:,'Load','ActivePower']   = LoadActivePower
Net.UBProfiles['Inputs',:,'Load','ReactivePower'] = LoadReactivePower
                                              
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
KKT[:J_active.shape[1],:J_active.shape[1]] = H+1e0*np.eye(J_active.shape[1])
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

dX = Net.VOptDispatch(dX)#/(Sol.cat+1e-12))

U,S,V = np.linalg.svd(J_active)
print "Min Singular value of J mat", np.min(S)
print "Max Singular value of J mat", np.max(S)
print "Conditioning of J (log scale)", np.log(np.max(S)/np.min(S))


V = V[:S.shape[0],:]
Sinv = np.diag(1/S)
S = np.diag(S)

np.allclose(J_active,np.dot(U,np.dot(S,V)))

np.dot(V.T,np.dot(Sinv,U.T))

#Net.ExtractInfo(dX, PlantPower = True, BusPower = True, TotalPower = True)
#Net.DYNSolvePlot(dX, dt = 1)

#E,W = np.linalg.eig(KKT)
#E = np.real(E)
#print "Min eigenvalue of KKT mat", np.min(np.abs(E))
#
#
#
#NJ = null(J_active.T,eps=1e-2)
##print "Null space: ", NJ

#SmallS = 10*np.min(S)

#ising = []
#for col in [0]:
#    for k in range(len(NJ[:,col])):
#        if (np.abs(NJ[k,col]) > 1e-1):
#            ising.append(k)
#        
#
#Jsing = J_active[ising,:]
#plt.subplot(2,1,1)
#plt.spy(Jsing)
#plt.subplot(2,1,2)
#plt.spy(J)
#plt.show()


Net.ExtractInfo(Sol, PlantPower = True, BusPower = True, TotalPower = True)
Net.DYNSolvePlot(Sol, dt = 1)
    
    
 
