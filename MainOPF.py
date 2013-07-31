# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien


"""

import os as os
import numpy as np
import sys
sys.path.append('/Users/sebastien/Desktop/DynOPFlow')
import DynOPFlow3
reload(DynOPFlow3)
from DynOPFlow3 import *

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

#Connection pair (i,j), Impedance
Graph = [[0,1,1+10j],
         [1,2,1+10j],
         [1,3,10+100j],
         [3,5,1+10j],
         [3,4,1+10j]] # Undirected connectivity graph

#Define Net properties

Net = PowerGrid(NBus,Graph)



Net.Flow(OPFSolver = True)



Net.PowerFlowBounds = {'Vmin' :           [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                       'Vmax' :           [10,  10,  10,   10,  10,  10],
                       'LineCurrentMax' : [inf, inf, inf, inf, inf] }

GridSetup = [{'Bus': 0, 'Property': 'PQ', 'P': -1., 'Q' : -0.75}]#,
             #{'Bus': 1, 'Property': 'Trans'},
             #{'Bus': 3, 'Property': 'Trans'}]


Sol = Net.OPFSolve(GridSetup)



Net._HessOPF.setInput(Net.OPF.output('x'),0)
Net._HessOPF.setInput(1.,1)
Net._HessOPF.setInput(1.,2)
Net._HessOPF.setInput(Net.OPF.output('lam_g'),3)
Net._HessOPF.evaluate()
H = Net._HessOPF.output()

Net._JacOPF.setInput(Net.OPF.output('x'),0)
Net._JacOPF.evaluate()
J = Net._JacOPF.output()

#g_sol = Net.gOPF(Net.OPF.output('g'))

g_sol = Net.OPF.output('g')
V_sol = Net.VOPF.cat
lam_g_sol = Net.OPF.output('lam_g')
i_active = []
g_active = []
for i in range(g_sol.shape[0]):
    if (Net.ubg.cat[i]-g_sol[i] < lam_g_sol[i]) or (g_sol[i]-Net.lbg.cat[i] < lam_g_sol[i]):
        i_active.append(i)
        print "Active Const: ", Net.gOPF.getLabel(i_active[-1])
        if (Net.ubg.cat[i]-g_sol[i] < lam_g_sol[i]):
            g_active.append(Net.ubg.cat[i]-g_sol[i])
        else:
            g_active.append(g_sol[i]-Net.lbg.cat[i])   


i_bound = []
Jbound = []
for i in range(Sol.shape[0]):
    if (Sol.cat[i] - Net.lbV.cat[i]  < 1e-10) or (Net.ubV.cat[i] - Sol.cat[i] < 1e-10):
        i_bound.append(i)
        print "Active Bound: ", Net.VOPF.getLabel(i_bound[-1])
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

Net._HessOPF.setInput(Net.OPF.output('x'),0)
#Net._HessOPF.setInput(Net.OPF.output('x'),0)
#Net._HessOPF.setInput(Net.OPF.output('x'),0)
Net._HessOPF.setInput(Net.OPF.output('lam_g'),3)
Net._HessOPF.evaluate()
H = Net._HessOPF.output()

Net._JacOPF.setInput(Net.OPF.output('x'),0)
Net._JacOPF.evaluate()
J = Net._JacOPF.output()


J_active = [np.array(J[i_active,:])]
J_active.append(Jbound)
J_active = np.concatenate(J_active,axis=0)

KKT = np.zeros([J_active.shape[1] + J_active.shape[0],J_active.shape[1] + J_active.shape[0]])
KKT[:J_active.shape[1],:J_active.shape[1]] = H+1e0*np.eye(J_active.shape[1])
KKT[J_active.shape[1]:,:J_active.shape[1]] = J_active
KKT[:J_active.shape[1],J_active.shape[1]:] = J_active.T

Net._JacCostOPF.setInput(Net.OPF.output('x'),0)
Net._JacCostOPF.evaluate()
JCost = Net._JacCostOPF.output()

RHS = np.zeros([J_active.shape[1] + J_active.shape[0],1])
RHS[:J_active.shape[1],:] = JCost
RHS[J_active.shape[1]:,:] = g_active

dX = np.linalg.solve(KKT,RHS)[:J_active.shape[1]]
print "Norm 2 of dX:", np.sqrt(mul(dX.T,dX))

dX = Net.VOPF(dX)#/(Sol.cat+1e-12))


E,W = np.linalg.eig(KKT)
E = np.real(E)
print "Min eigenvalue of KKT mat", np.min(np.abs(E))

U,S,V = np.linalg.svd(J_active)
print "Min Singular value of J mat", np.min(S)



NJ = null(J_active.T,eps=1e-2)
#print "Null space: ", NJ

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



    
    
 


