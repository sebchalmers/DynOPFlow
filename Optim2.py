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

#Define Net properties
Net = PowerGrid(NBus,Graph)
Net.Flow()

Net.PowerFlowBounds = {
                       'Vmin' :           [0.25 for k in range(NBus)],
                       'Vmax' :           [10.0 for k in range(NBus)],
                       'LineCurrentMax' : [inf for k in range(5)    ]
                      }


dt = 1.



#####   Thermal plants  #####
plantKeys = []
for bus in [1,2,3,4,5]:
    plantKeys.append('Thermal'+str(bus))
    Thermal           = Plant(Bus = bus, R = 0.0, label = plantKeys[-1])
    ThermalPower      = Thermal.Inputs['Power']
    ThermalPower_prev = Thermal.InputsPrev['Power']

    Cost =  (ThermalPower - ThermalPower_prev)**2 + 1e3*ThermalPower 
    Thermal.setCost(Cost)

    Net.addPlant(Thermal)

    Thermal.UB['Inputs','Power'] = 1.           


#####   Load      ######
Load = Plant(Load = True, Bus = 0, label = 'Load')
Net.addPlant(Load)

Load.LB['Inputs',  'ActivePower'] = -1.
Load.LB['Inputs','ReactivePower'] =  -0.75
Load.UB['Inputs',  'ActivePower'] = -1.
Load.UB['Inputs','ReactivePower'] =  -0.75

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
        
Horizon    = 48
#Nsimulation = int(10*24)

Net.Profiles(Horizon)

Nprofile = Net.Nprofile

LoadActivePower   = [0.3*np.cos(2*np.pi*k*dt/24.) - 1.  for k in range(Nprofile)]
LoadReactivePower = [0.75*LoadActivePower[k] for k in range(Nprofile)]

Net.Dispatch(Horizon = Horizon),# Simulation = Nsimulation)

#Initial conditions (set inf in x0 to free the initial conditions)
u0 = Net.u0()

for plant in plantKeys:
    u0[plant,'Power']  = 0.


#Make initial guess
init = Net.init()

                
Net.LBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.LBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
Net.UBInputProfiles['Load',:,'ActivePower']   = LoadActivePower
Net.UBInputProfiles['Load',:,'ReactivePower'] = LoadReactivePower
                                              
Sol,_ = Net.DYNSolve(u0 = u0, init = init)

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
        print "Active Const: ", Net.gOptDispatch.getLabel(i_active[-1])
        if (Net.ubg.cat[i]-g_sol[i] < lam_g_sol[i]):
            g_active.append(Net.ubg.cat[i]-g_sol[i])
        else:
            g_active.append(g_sol[i]-Net.lbg.cat[i])   


i_bound = []
Jbound = []
for i in range(Sol.shape[0]):
    if (Sol.cat[i] - Net.lbV.cat[i]  < 1e-10) or (Net.ubV.cat[i] - Sol.cat[i] < 1e-10):
        i_bound.append(i)
        print "Active Bound: ", Net.VOptDispatch.getLabel(i_bound[-1])
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

for i in range(Sol.shape[0]):
    print Sol.getLabel(i),"=", Sol.cat[i]
    
dX = np.linalg.solve(KKT,RHS)[:J_active.shape[1]]
print "Norm 2 of dX:", np.sqrt(mul(dX.T,dX))

dX = Net.VOptDispatch(dX)#/(Sol.cat+1e-12))




##Net.ExtractInfo(dX, PlantPower = True, BusPower = True, TotalPower = True)
##Net.DYNSolvePlot(dX, dt = 1)
#
#E,W = np.linalg.eig(KKT)
#E = np.real(E)
#print "Min eigenvalue of KKT mat", np.min(np.abs(E))
#
U,S,V = np.linalg.svd(J_active)
print "Min Singular value of J mat", np.min(S)
print "Max Singular value of J mat", np.max(S)
print "Conditioning of J (log scale)", np.log(np.max(S)/np.min(S))
#
#NJ = null(J_active.T,eps=1e-2)
##print "Null space: ", NJ

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



    
    
 
