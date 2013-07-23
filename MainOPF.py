# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien


"""

import DynaPFlow
reload(DynaPFlow)

from DynaPFlow import *


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



Net.Flow(OPFSolver = 'True')



Net.PowerFlowBounds = {'Vmin' :           [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                       'Vmax' :           [10,  10,  10,   10,  10,  10],
                       'LineCurrentMax' : [inf, inf, inf, inf, inf] }

GridSetup = [{'Bus': 0, 'Property': 'slack', 'V': 1.},
             {'Bus': 2, 'Property': 'PQ', 'P': 1., 'Q': 0.}]


OPF1 = Net.OPFSolve(GridSetup)



Net.PowerFlowBounds = {}
GridSetup = [{'Bus': 0, 'Property': 'slack', 'V': 1.},
             {'Bus': 2, 'Property': 'PQ', 'P': 1., 'Q': 0.8},        
             {'Bus': 4, 'Property': 'PV', 'P': 1., 'V': 1.2}]

OPF2 = Net.OPFSolve(GridSetup)

