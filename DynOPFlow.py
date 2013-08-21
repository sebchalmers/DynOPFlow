# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author:

Sebastien Gros
Assistant Professor
 
Department of Signals and Systems
Chalmers University of Technology
SE-412 96 Göteborg, SWEDEN
grosse@chalmers.se

Python/casADi Module:
NMPC for Dynamic Optimal Power Flow and Power Dispatch

Requires the installation of the open-source Python module casADi together with the NLP solver ipopt

Required version of CasADi: v1.7.x

"""

from casadi import *
from casadi.tools import *
import math as math

import numpy as np
import matplotlib.pyplot as plt
import random as rand
#from   ODE import *

plt.rcParams['text.usetex'] = False
    
#Fixed constants
#rho_air   = 1.2
#rho_water = 1e3
#gravity   = 9.81
    
def assertList(var):
    if not(isinstance(var,list)):
        var = [var]
    return var        

class Plant:
    def __init__(self, Inputs = [], States = [], ExtParameters = [], R = 0., Directionality = 'Mono', Load = False, Bus = [], label = []):        # set of reasons (strings) why the Dae cannot be modified (add new x/z/u/p/output)

        
        self._frozen = False
        self.Bus     = Bus
        self.label = label
        self.Directionality = Directionality
        self._Load = Load
        self.R = R
        
        #Plant Default Input structure (list ready to be embedded in a struct_sym)        
        InputList = [entry("CurrentReal"),
                     entry("CurrentImag")]
        
        #Structure for INPUTS of various power plants
        if      (Load == True):
            InputList.append(entry("ActivePower"))
            InputList.append(entry("ReactivePower"))
        elif    (Directionality == 'Bi'):
            InputList.append(entry("Pcharge"))
            InputList.append(entry("Pdischarge"))
        elif    (Directionality == 'Mono'):
            InputList.append(entry("Power"))
        else:
            print "Illegal option, ignored"
            return
    

        if (len(Inputs) > 0):
            self._additionalInputs = Inputs #Keep the list for plotting purposes
            for key in assertList(Inputs):
                InputList.append(entry(key))
    
        # States declared by the user
        #if (len(States) > 0): ################ INTRODUCE THIS ######## !!!!
        StateList = []
        for key in assertList(States):
            StateList.append(entry(key)) 
    
        # External parameters declared by the user
        if (len(ExtParameters) > 0):
            ExtParamList = []
            for key in assertList(ExtParameters):
                ExtParamList.append(entry(key))
            self.ExtParameters      = struct_ssym(ExtParamList)    

        # lists of names (strings)
        self.States            = struct_ssym(StateList)
        self.Inputs            = struct_ssym(InputList)
        self.InputsPrev        = struct_ssym(InputList)
        
        
        #Structure for plant bounds
        Bound = [entry('Inputs',   struct = self.Inputs)]
        if (len(self.States.keys()) > 0):
            Bound.append(entry('States',   struct = self.States))            
        Bound = struct_ssym(Bound)
        
        self.LB = Bound(-inf)
        self.UB = Bound( inf)
        
        if (Directionality == 'Mono') and (Load == False):
            self.LB['Inputs','Power'] = 0.
        elif (Directionality == 'Bi') and (Load == False):
            self.LB['Inputs','Pcharge'] = 0.
            self.LB['Inputs','Pdischarge'] = 0.
        else:
            self.LB['Inputs','ActivePower'] = 0.
            

            
    def setDynamics(self, RHS = [], dt = 1., nstep = 10):
        if (self._frozen == True):
            print "Plant already added to the grid, call ignored"
            return

        
        print "Right-Hand side: ", RHS
        
        if isinstance(RHS,list):
            RHS = veccat(RHS)
        
        X = self.States
        U = self.Inputs
        
        dtRK4 = dt/float(nstep)
        
        fimplicit = SXFunction(daeIn(x=X,p=U),daeOut(ode=RHS))
        fimplicit.init()
        
        [k1]  = daeOut(fimplicit.eval(daeIn(  x=X,p=U                )),"ode")
        [k2]  = daeOut(fimplicit.eval(daeIn(  x=X+0.5*dtRK4*k1,p=U   )),"ode")
        [k3]  = daeOut(fimplicit.eval(daeIn(  x=X+0.5*dtRK4*k2,p=U   )),"ode")
        [k4]  = daeOut(fimplicit.eval(daeIn(  x=X+dtRK4*k3,p=U       )),"ode")
        
        rk4_step = SXFunction([X,U],[X + (1./6)*dtRK4*(k1 + 2*k2 + 2*k3 + k4)])
        rk4_step.init()
        
        #CONSTRUCT SHOOTING
        # -----------------------------------------------------------------------------
        out = X
        for i in range(0,nstep):
            [out] = rk4_step.eval([out,U])
            
        Shoot = SXFunction([X,U],[out])        
           
        #Shoot = SXFunction([X,U],[X+dt*RHS])  #Invoke a 1st order Euler... 
        
        Shoot.init()
    
        self._Shoot = Shoot
     
    def _BuildFunc(self, Expr, Terminal):
    
        X      = self.States
        U      = self.Inputs
        Uprev  = self.InputsPrev
                   
          
        if Terminal == False:
            listFuncInput = [U, Uprev]            
            if    (X.size > 0):
                listFuncInput.append(X)
        else:
            listFuncInput = [X]

        if hasattr(self,'ExtParameters'):
            listFuncInput.append(self.ExtParameters)
            
        Func = SXFunction(listFuncInput,[Expr])
        Func.init()
        
        return Func
    
    def setConstraints(self, Const, Terminal = False):
        
        if (self._frozen == True):
            print "Plant already added to the grid, call ignored"
            return
        
        if not(isinstance(Const,list)):
            Const = [Const]
        
        #ConstFunc = self._BuildFunc(veccat(Const), Terminal)
                
        if    (Terminal == False):
            self._StageConst    = self._BuildFunc(veccat(Const), Terminal)
        elif  (Terminal == True): 
            self._TerminalConst = self._BuildFunc(veccat(Const), Terminal)        
                
        
    def setCost(self, Cost, Terminal = False):
        if (self._frozen == True):
            print "Plant already added to the grid, call ignored"
            return

        #CostFunc = self._BuildFunc(Cost, Terminal)
        
        if    (Terminal == False):
            self._StageCost    = self._BuildFunc(Cost, Terminal)
        elif  (Terminal == True): 
            self._TerminalCost = self._BuildFunc(Cost, Terminal)
     


class PowerGrid:
    """
    Generates:
    - Power FLow equations and Power plants dynamics
    - OPF solver
    - optimal grid control solver
    """ 

    def __init__(self, NBus = 0,Graph = []):        
        self.NBus = NBus
        self.Graph = Graph
        self.PlantList = []
        
        self._hasStates = False
        
        self.PowerFlowBounds = {'Vmin' :           0,
                                'Vmax' :           inf,
                                'LineCurrentMax' : inf }
    

    #CONSTRUCT POWER FLOW    
    def Flow(self, OPFSolver = False):
        NBus =  self.NBus
        NLine = np.size(    self.Graph       ,axis = 0)

        print "Constructing Power Flow Equations, #Bus =",NBus, ", #Line =",NLine
       
        #CONSTRUCT THE POWER FLOW EQUATIONS

        Graph = self.Graph
        
        # Bus admittance matrix: Inodal_injection = Y*V (Nodal current injection)
        Y = np.array([ [ 0.0*1j for i in range(NBus) ] for j in range(NBus) ])
        for k in range(NLine):
            Y[Graph[k][0],Graph[k][0]] += 1/Graph[k][2]
            Y[Graph[k][1],Graph[k][1]] += 1/Graph[k][2]
            Y[Graph[k][0],Graph[k][1]] -= 1/Graph[k][2]
            Y[Graph[k][1],Graph[k][0]] -= 1/Graph[k][2]
            
        # Line admittance matrix (directed): Iline = L*V  
        L = np.array([ [ 0.0*1j for i in range(NBus) ] for j in range(NLine) ])
        for k in range(NLine):
            L[k,Graph[k][0]] =  1/Graph[k][2]
            L[k,Graph[k][1]] = -1/Graph[k][2]
           
        ######## BUILD POWER FLOW EQUATIONS (Results in Function: Bus voltage, Bus power -> Residual, to be satisfied at every time stage)
        
        #Bus Voltages (real and complex parts)
        BusVoltages = struct_ssym([entry("Real",repeat = NBus),
                                   entry("Imag",repeat = NBus)])
                                            
        #Bus currents  (for current limitation) I = Y*V
        BusCurrentsReal = mul(np.real(Y),BusVoltages["Real",veccat]) - mul(np.imag(Y),BusVoltages["Imag",veccat]) 
        BusCurrentsImag = mul(np.real(Y),BusVoltages["Imag",veccat]) + mul(np.imag(Y),BusVoltages["Real",veccat]) 
        BusCurrents2 = BusCurrentsReal*BusCurrentsReal + BusCurrentsImag*BusCurrentsImag
        
        #Bus voltage modules square (for voltage limitation)
        BusVoltages2 = BusVoltages["Real",veccat]*BusVoltages["Real",veccat] + BusVoltages["Imag",veccat]*BusVoltages["Imag",veccat]
        
        #Line currents
        LineCurrentsReal = mul(np.real(L),BusVoltages["Real",veccat]) - mul(np.imag(L),BusVoltages["Imag",veccat]) 
        LineCurrentsImag = mul(np.real(L),BusVoltages["Imag",veccat]) + mul(np.imag(L),BusVoltages["Real",veccat]) 
        LineCurrents2 = LineCurrentsReal*LineCurrentsReal + LineCurrentsImag*LineCurrentsImag
        
        
        #Build complex power injections at the bus
        SReal = BusCurrentsReal*BusVoltages["Real",veccat] + BusCurrentsImag*BusVoltages["Imag",veccat]
        SImag = BusCurrentsReal*BusVoltages["Imag",veccat] - BusCurrentsImag*BusVoltages["Real",veccat] 
        
        #Create functions for Current and Voltages**2
        self.BusActivePowerFunc   = SXFunction([BusVoltages],[SReal])
        self.BusReactivePowerFunc = SXFunction([BusVoltages],[SImag])
        self.BusCurrentsRealFunc = SXFunction([BusVoltages],[BusCurrentsReal])
        self.BusCurrentsImagFunc = SXFunction([BusVoltages],[BusCurrentsImag])
        self.LineCurrents2Func   = SXFunction([BusVoltages],[LineCurrents2])
        self.BusVoltages2Func   = SXFunction([BusVoltages],[BusVoltages2])
        
        self.BusVoltages = BusVoltages
        self.BusActivePowerFunc.init()
        self.BusReactivePowerFunc.init()
        self.BusCurrentsRealFunc.init()
        self.BusCurrentsImagFunc.init()
        self.LineCurrents2Func.init()
        self.BusVoltages2Func.init()
        self.OPF = True

        #CONSTRUCT OPF SOLVER IS ASKED IN THE FLOW FUNCTION OPTIONS 
        if (OPFSolver == True):
            print "Construct OPF Solver"
            
            Power = struct_ssym([
                                    entry('Active',   repeat = NBus),
                                    entry('Reactive', repeat = NBus)
                                ])
            
            V = struct_msym([
                                    entry('BusPower', struct = Power),
                                    entry('BusVoltages', struct = BusVoltages)
                            ])
            
            [BusActivePower] =   self.BusActivePowerFunc.call([V['BusVoltages']])
            [BusReactivePower] = self.BusReactivePowerFunc.call([V['BusVoltages']])
            [LineCurrents2] =    self.LineCurrents2Func.call([V['BusVoltages']])
            [BusVoltages2] =     self.BusVoltages2Func.call([V['BusVoltages']])
            
            ActivePowerBalance =     BusActivePower - V['BusPower','Active',  veccat]
            ReactivePowerBalance = BusReactivePower - V['BusPower','Reactive',veccat]

            g = struct_MX([
                              entry('ActivePower',   expr = ActivePowerBalance),
                              entry('ReactivePower', expr = ReactivePowerBalance),
                              entry('LineCurrents2', expr = LineCurrents2),
                              entry('BusVoltages2',  expr = BusVoltages2)
                          ])
            
            Cost = 0
            for line in range(NLine):
                Cost += np.real(Graph[line][2])*LineCurrents2[line]
            
            nl = MXFunction(nlpIn(x=V),nlpOut(f=Cost,g=g))
                   
            nl.init()
            
            # set-up solver
            solver = IpoptSolver(nl)
            #solver.setOption("print_level",0)
            solver.setOption("expand",True)
            solver.setOption("parametric",False)    
            solver.setOption("generate_hessian",True)
            solver.setOption("max_iter",1000)
            solver.setOption("tol",1e-6)
            solver.setOption("linear_solver","ma27")
            
            solver.init()
            Hessian = solver.hessLag()
            Hessian.init()
            
            Jacobian = solver.jacG()
            Jacobian.init()
    
            JacCost = solver.gradF()
            JacCost.init()
        
            self.VOPF     = V
            self.gOPF = g
            self.OPF = solver
                 
            self._HessOPF = Hessian
            self._JacOPF = Jacobian
            self._JacCostOPF  = JacCost
            
    def OPFSolve(self, Grid = []):
        
        lbV  =  self.VOPF(-inf)
        ubV  =  self.VOPF( inf)
        lbg  =  self.gOPF()
        ubg  =  self.gOPF()
        
        if not(hasattr(self,'OPF')):
            #Check that .Flow() has been called
            print "You must call .Flow(OPFSolver = True) to setup OPF before calling .OPFSolve()"
            return []

        if (self.OPF == True):
            #Check that a solver exists
            print "You must call .Flow(OPFSolver = True) to setup OPF before calling .OPFSolve()"
            return []
    
    
        #Set initial guess
        init =  self.VOPF()
        init['BusVoltages','Real',veccat] = 1.
        
        #Set the bounds (default values if not defined)
        lbV  =  self.VOPF(-inf)
        ubV  =  self.VOPF( inf)
        lbg  =  self.gOPF()
        ubg  =  self.gOPF()

        #Ascertain the completness of the PowerFlowBounds dictionary, complete if necessary
        if not('Vmin' in self.PowerFlowBounds.keys()):
            self.PowerFlowBounds['Vmin'] = 0
            print "Min Bus Voltage not provided, default value assigned (0)"
        if not('Vmax' in self.PowerFlowBounds.keys()):
            self.PowerFlowBounds['Vmax'] = inf
            print "Max Bus Voltage not provided, default value assigned (inf)"
        if not('LineCurrentMax' in self.PowerFlowBounds.keys()):
            self.PowerFlowBounds['LineCurrentMax'] = inf
            print "Max Line Current not provided, default value assigned (inf)"


        ubg['LineCurrents2'] = np.array(self.PowerFlowBounds['LineCurrentMax'])**2
        lbg['BusVoltages2']  = np.array(self.PowerFlowBounds['Vmin'])**2
        ubg['BusVoltages2']  = np.array(self.PowerFlowBounds['Vmax'])**2




        #Assign Network operational conditions
        for entry in range(np.size(Grid)):                
            if      (Grid[entry]['Property'] == 'slack'):
                print "Bus", Grid[entry]['Bus'],"is slack"
                lbg['BusVoltages2',Grid[entry]['Bus']]       = Grid[entry]['V']**2
                ubg['BusVoltages2',Grid[entry]['Bus']]       = Grid[entry]['V']**2
                lbV['BusVoltages','Imag',Grid[entry]['Bus']] =       0.0
                ubV['BusVoltages','Imag',Grid[entry]['Bus']] =       0.0
                
            elif    (Grid[entry]['Property'] == 'PV'):
                print "Bus", Grid[entry]['Bus'],"is PV"
                lbg['BusVoltages2',Grid[entry]['Bus']]       = Grid[entry]['V']**2
                ubg['BusVoltages2',Grid[entry]['Bus']]       = Grid[entry]['V']**2
                lbV['BusPower','Active',Grid[entry]['Bus']]  = Grid[entry]['P']
                ubV['BusPower','Active',Grid[entry]['Bus']]  = Grid[entry]['P']

            elif    (Grid[entry]['Property'] == 'PQ'):
                print "Bus", Grid[entry]['Bus'],"is PQ"
                lbV['BusPower','Active',Grid[entry]['Bus']]    = Grid[entry]['P']
                ubV['BusPower','Active',Grid[entry]['Bus']]    = Grid[entry]['P']
                lbV['BusPower','Reactive',Grid[entry]['Bus']]  = Grid[entry]['Q']
                ubV['BusPower','Reactive',Grid[entry]['Bus']]  = Grid[entry]['Q']
            
        
        
        self.OPF.setInput( lbV,     "lbx")
        self.OPF.setInput( ubV,     "ubx")
        self.OPF.setInput(init,     "x0" )
        self.OPF.setInput( lbg,     "lbg")
        self.OPF.setInput( ubg,     "ubg")
        
        self.OPF.solve()
         
        self.lbg = lbg
        self.ubg = ubg
        self.lbV = lbV
        self.ubV = ubV
        
        return  self.VOPF(self.OPF.output('x'))

###########     POWER DISPACTH PROBLEM   ##########

    def addPlant(self, plant):
        if isinstance(plant,list): #Treat list of plants
            for plant_k in plant:
                self.addPlant(plant_k)
        else:        
            if (plant._frozen == True):
                print "Plant already added to the grid, call ignored"
                return
            
            self.PlantList.append(plant)
            plant._frozen = True
            if hasattr(plant,'_Shoot'):
                self._hasStates = True
                
    def _VariableConstructor(self, N):
        ###### CONSTRUCT DECISION VARIABLES OF LENGTH N #######
        List = []
        for plant in self.PlantList:
            List.append(entry(plant.label,  struct = plant.Inputs))
        Inputs = struct_ssym(List)
        
        List = []
        for plant in self.PlantList:
            if (len(plant.States.keys()) > 0):
                List.append(entry(plant.label,  struct = plant.States))
        States = struct_ssym(List)
        
        
        #Structures to manipulate initial conditions and inputs
        u0 = struct_msym(Inputs)
        x0 = struct_msym(States)
        
        #User-specified additional parameters
        EPList = []
        for plant in self.PlantList:
            if hasattr(plant,'ExtParameters'):
                EPList.append(entry(plant.label, struct = plant.ExtParameters))
        
        ExtParameters = struct_msym(EPList)
        
        EP = struct_msym([
                            entry('u0',              struct = u0),
                            entry('ExtParameters',   struct = ExtParameters)
                         ])
        
        Vlist = []
        Vlist.append(entry("BusVoltages",             repeat = N,     struct = self.BusVoltages))
        if (self._hasStates == True):
            Vlist.append(entry("States",              repeat = N+1,   struct = States))
        Vlist.append(entry("Inputs",                  repeat = N,     struct = Inputs))
        
        V = struct_msym(Vlist)

        return V, u0, x0, EP, ExtParameters
    

    def _CostConstructor(self, V, EP, Nstage, GridLoss):
        """
        Constructor for the Cost function, handy to build the cost for different V:s
        """
        
        Cost_Lagrange = 0
        #Grid loss
        if (GridLoss == True):
            NLine     =  len(    self.Graph   )        
            for k in range(Nstage):
                # Grid losses
                [LineCurrents2_k] = self.LineCurrents2Func.call([V['BusVoltages',k]])
                for line in range(NLine):
                    Cost_Lagrange += np.real(self.Graph[line][2])*LineCurrents2_k[line]
        
        #Plants Lagrange cost
        for plant in self.PlantList:
            for k in range(Nstage):    
                if (hasattr(plant,'_StageCost')):
                    CostInputList = [V['Inputs',k,plant.label]]
                    if (k==0):
                        CostInputList.append(  EP['u0',plant.label]       )
                    else:
                        CostInputList.append(  V['Inputs',k-1,plant.label])

                    if  (plant.States.size > 0):
                        CostInputList.append(V['States',k,plant.label])

                    if hasattr(plant,'ExtParameters'):
                        CostInputList.append(EP['ExtParameters',plant.label])

                    [Cost_k] = plant._StageCost.call(CostInputList)
                    Cost_Lagrange += Cost_k
        
        #Plants Terminal cost
        Cost_Terminal = 0
        for plant in self.PlantList:            
            if (hasattr(plant,'_TerminalCost')):
                CostInputList = [V['States',-1,plant.label]]
                if hasattr(plant,'ExtParameters'):
                        CostInputList.append(EP['ExtParameters',plant.label])
                        
                [Cost_k] = plant._TerminalCost.call(CostInputList)
                Cost_Terminal += Cost_k
                        
        Cost = (Cost_Lagrange+Cost_Terminal)/Nstage
        
        LagrangeCostFunc = MXFunction([V,EP],[Cost_Lagrange])
        LagrangeCostFunc.init() 
        TerminalCostFunc = MXFunction([V,EP],[Cost_Terminal])
        TerminalCostFunc.init()
        
        return Cost, LagrangeCostFunc, TerminalCostFunc

    def Dispatch(self, Horizon = 24, Simulation = 0, GridLoss = True):
        """
        Constructs the power dispatch problem, default Horizon length (if argument Horizon is not provided) is 24 time units
        """
        if (self.OPF == False):
            Power.Flow(self)
            
        print "Construct Dynamic OPF"
        
        #THorizon = self.THorizon
        Nstage = Horizon#self.TimeSetup['Horizon'] 

        NBus = self.NBus
        NLine     =  len(    self.Graph   )


        TransferBus = []
        BusProperties = []
        for Bus in range(NBus):
            Busk = 'transfer'
            BuskProperties = []
            for plant in self.PlantList:  
                if ( plant.Bus == Bus):
                    Busk = 'open'
                    BuskProperties.append(plant.label)
            if (BuskProperties == []):
                TransferBus.append(Bus)               
            BusProperties.append({Bus: BuskProperties})           
            
  
        ###################   CONSTRUCT VARIABLES    ########################
             
        V, u0, x0, EP, ExtParameters = self._VariableConstructor(Nstage)
        
        
        #Structure for storing NMPC solutions if Nsim provided
        if (Simulation > 0):
             
            Vstore,_,_,_,_ = self._VariableConstructor(Simulation)
            self.Vstore = Vstore()
                  
        ###############################     BUILD COST AND CONSTRAINTS         ###############################
        
        Cost, LagrangeCostFunc, TerminalCostFunc = self._CostConstructor(V, EP, Nstage, GridLoss)
        
        if (Simulation > 0):   
            _, self.LagrangeCost, self.TerminalCost = self._CostConstructor(Vstore, EP, Simulation, GridLoss)
        else:
            self.LagrangeCost = LagrangeCostFunc
            self.TerminalCost = TerminalCostFunc

        # OPF constraints
        CurrentBalance  = []
        LineCurrents2   = []
        BusVoltages2    = []
        
        # Generic constraints
        PeriodicConst   = [] 
        EquConst        = []
        IneqConst       = []
        
        # Thermal constraints
        ThermalConst    = []
        ThermalConstExt = []
       

        
    
        #########    BUILD COST & CONSTRAINTS    #######
        for k in range(Nstage): #k is reserved for time instant throughout the code
            
            ### CONSTRUCT POWER FLOW      

             #Construct (Bus Voltages)**2 and (Line current)**2 for bounding module 
            [LineCurrents2_k] = self.LineCurrents2Func.call([V['BusVoltages',k]])
            [BusVoltages2_k]  =  self.BusVoltages2Func.call([V['BusVoltages',k]])
            LineCurrents2.append(LineCurrents2_k)
            BusVoltages2.append(BusVoltages2_k)
            
            #Compute Bus Injection Currents
            [CurrentsBalanceReal] = self.BusCurrentsRealFunc.call([V['BusVoltages',k]])
            [CurrentsBalanceImag] = self.BusCurrentsImagFunc.call([V['BusVoltages',k]])       
        
            
            for plant in self.PlantList:

                #Bus Voltage for the selected plant/load
                BusVoltageReal = V['BusVoltages',k,'Real'][plant.Bus]
                BusVoltageImag = V['BusVoltages',k,'Imag'][plant.Bus]
    
                #Plant Current of the selected plant/load
                PlantCurrentReal = V['Inputs',k,plant.label,'CurrentReal']
                PlantCurrentImag = V['Inputs',k,plant.label,'CurrentImag']
                
                # Balance the participating currents of the various plants and loads with
                # the current injection @ the corresponding buses

                CurrentsBalanceReal[plant.Bus] -= PlantCurrentReal
                CurrentsBalanceImag[plant.Bus] -= PlantCurrentImag
                           
                # Re{V.iplant*} -> "Participating Active Power" // Im{V.iplant*} -> "Participating Reactive Power"
                ParticipatingActivePower   = BusVoltageReal*PlantCurrentReal + BusVoltageImag*PlantCurrentImag
                ParticipatingReactivePower = BusVoltageImag*PlantCurrentReal - BusVoltageReal*PlantCurrentImag
                
                # Plant participating current squared, i.e. |i|**2
                PlantCurrent2 = PlantCurrentReal*PlantCurrentReal + PlantCurrentImag*PlantCurrentImag
                
                if (plant._Load == True):
                    #Load fixing: [Active, Reactive] = Consumed Active / Reactive power
                    EquConst.append(ParticipatingActivePower   - V['Inputs',k,plant.label,'ActivePower'])
                    EquConst.append(ParticipatingReactivePower - V['Inputs',k,plant.label,'ReactivePower'])
                else:
                    if (plant.Directionality == 'Mono'):
                        PlantPower = V['Inputs',k,plant.label,'Power']
                    else:
                        PlantPower = V['Inputs',k,plant.label,'Pdischarge'] - V['Inputs',k,plant.label,'Pcharge']  
        
                    #Compute balance between Pmech and participating power for each plant
                    # ParticipatingPower + R*|iplant|**2 - PlantPower = 0
                    EquConst.append(ParticipatingActivePower + plant.R*PlantCurrent2 - PlantPower)
        
            CurrentBalance.append(CurrentsBalanceReal)
            CurrentBalance.append(CurrentsBalanceImag)
        
            ### CONSTRUCT DYNAMIC CONSTRAINTS
            
            for plant in self.PlantList:
                if hasattr(plant,'_Shoot'):
                    [Xp] = plant._Shoot.call([V['States',k,plant.label],V['Inputs',k,plant.label]])    
                    EquConst.append(Xp-V['States',k+1,plant.label])
            
                #A bit ugly...
                if hasattr(plant,'_StageConst'):
                    #print "Plant", plant.label, "has stage inequality constraints"
                    ConstInputList = [V['Inputs',k,plant.label]]
                    if (k==0):
                        ConstInputList.append(  EP['u0',plant.label]       )
                    else:
                        ConstInputList.append(  V['Inputs',k-1,plant.label])
            
                    if  (plant.States.size > 0):
                        ConstInputList.append(V['States',k,plant.label])
            
                    [Const_k] = plant._StageConst.call(ConstInputList)
                    IneqConst.append(Const_k)
    

        ### END OF STAGE CONSTRAINTS
        for plant in self.PlantList:
            if (hasattr(plant,'_TerminalConst')):
                #print "Plant", plant.label, "has terminal inequality constraints"
                [Const_k] = plant._TerminalConst.call([V['States',-1,plant.label]])
                IneqConst.append(Const_k)

        
        ####  PERIODIC CONSTRAINTS
        #PeriodicConst.append(V['States',-1,'Storage','Energy']    - V['States',0,'Storage','Energy'])
        #PeriodicConst.append(V['States',-1,'Hydro','WaterHeight'] - V['States',0,'Hydro','WaterHeight'])
        #PeriodicConst.append(V['Inputs',-1,'Thermal','Power']     - V['Inputs',0,'Thermal','Power'])

        
        ######## END CONSTRAINTS BUILDING ######
        


        g = struct_MX([
          entry("CurrentBalance", expr = CurrentBalance),
          entry("BusVoltages2",   expr = BusVoltages2),
          entry("LineCurrents2",  expr = LineCurrents2),
          #entry('Periodic',       expr = PeriodicConst),
          entry('EquConst',       expr = veccat(EquConst)),
          entry('IneqConst',      expr = veccat(IneqConst))
        ])
        
        nl = MXFunction(nlpIn(x=V,p=EP),nlpOut(f=Cost,g=g))
        nl.init()
        
        # set-up solver
        solver = IpoptSolver(nl)  
        solver.setOption("expand",True)
        solver.setOption("print_level",0)
        solver.setOption("parametric",True)  
        solver.setOption("hessian_approximation","exact")
        solver.setOption("max_iter",2000)
        solver.setOption("tol",1e-6)
        solver.setOption("linear_solver","ma27")
        
        solver.init()
        Hessian = solver.hessLag()
        Hessian.init()
        
        Jacobian = solver.jacG()
        Jacobian.init()
        
        JacCost = solver.gradF()
        JacCost.init()
        
        self._HessOptDispatch = Hessian
        self._JacOptDispatch  = Jacobian
        self._JacCostOptDispatch  = JacCost
                
        self.u0            = u0
        self.x0            = x0
        self.ExtParameters = ExtParameters
        self._EP           = EP
        self.VOptDispatch  = V
        self.OptDispatch   = solver        
        self.gOptDispatch  = g
        
        self.Properties = BusProperties

        print self.Properties
        ##############  SOLVER CONSTRUCTED  ##############
        
       
        
    
    #BUILD FIRST INITIAL GUESS   
    def init(self, x0 = [], u0 = []):
        
        init = self.VOptDispatch()
        
        NBus = self.NBus
        NLine     =  len(    self.Graph   )

        for plant in self.PlantList:
            if hasattr(plant,'_Shoot'):
                init['States',:,plant.label] = 0.5*(plant.LB['States'] + plant.UB['States'])

        for index in range(init.size):
            if not(init.cat[index] < inf):
                init.cat[index] = 0.
            
            
        for bus in range(NBus):
            init['BusVoltages',:,'Real',bus] = 0.5*(self.PowerFlowBounds['Vmin'][bus]+self.PowerFlowBounds['Vmax'][bus])
            init['BusVoltages',:,'Imag',bus] = 0.0
            
        init['Inputs',:,...,'CurrentReal'] = 1.0
        init['Inputs',:,...,'CurrentImag'] = 1.0   
        
        return init
        
    def Profiles(self, N = 0):
        """
        CREATE A STRUCTURE FOR HANDLING THE PROFILES OF THE POWER GRID:
        If no argument passed, the profiles have the horizon length, if argument N is assigned, the profiles have the length of NSample
        """
         
        if (N == 0):   
            if hasattr(self,'VOptDispatch'):
                Nstage   = len(self.VOptDispatch['Inputs'])
                Nprofile = Nstage + 1
            else:
                print "Profile Error: cannot resolve the length of profile. Specify a horizon (N = ...) or call .Dispatch first" 
                return
        else:
            if hasattr(self,'VOptDispatch'):
                Nstage   = len(self.VOptDispatch['Inputs'])
                Nprofile = N + Nstage + 1
            else:
                Nprofile = N + 1
                
        self.Nprofile = Nprofile
    
        VProfile,_,_,_,_ = self._VariableConstructor(self.Nprofile)  
    
        self.LBProfiles = VProfile()
        self.UBProfiles = VProfile()
        

        for plant in self.PlantList:
            self.LBProfiles['Inputs',:,plant.label] = plant.LB['Inputs']
            self.UBProfiles['Inputs',:,plant.label] = plant.UB['Inputs']
            if hasattr(plant,'_Shoot'):
                self.LBProfiles['States',:,plant.label] = plant.LB['States']
                self.UBProfiles['States',:,plant.label] = plant.UB['States']
    
        return Nprofile
    #ASSIGN PROFILES & SOLVE

    def DYNSolve(self, x0 = [], u0 = 0., ExtParameters = [], init = [], time = 0, Periodic = False):
        
        lbV  =  self.VOptDispatch(-inf)
        ubV  =  self.VOptDispatch( inf)
        lbg  =  self.gOptDispatch()
        ubg  =  self.gOptDispatch()
        
        ubg["IneqConst"] = 0.
        lbg["IneqConst"] = -inf

        NBus = self.NBus
        NLine     =  len(    self.Graph   )

        ####### SETUP THE BOUNDS #########            
        for plant in self.PlantList:
            lbV['Inputs',:,plant.label] = self.LBProfiles['Inputs',time:,plant.label]
            ubV['Inputs',:,plant.label] = self.UBProfiles['Inputs',time:,plant.label]
            if hasattr(plant,'_Shoot'):
                lbV['States',:,plant.label] = self.LBProfiles['States',time:,plant.label]
                ubV['States',:,plant.label] = self.UBProfiles['States',time:,plant.label]
                
        #Power flow limitations
        lbg["BusVoltages2"]  = np.array(self.PowerFlowBounds['Vmin'])**2
        ubg["BusVoltages2"]  = np.array(self.PowerFlowBounds['Vmax'])**2
        ubg["LineCurrents2"] = np.array(self.PowerFlowBounds['LineCurrentMax'])**2
        
        #Introduce additional bounds on all current and voltages (taken from Power flow limitation)
        # Bus voltages
        for bus in range(NBus):
            ubV['BusVoltages',:,'Real',bus] =  self.PowerFlowBounds['Vmax'][bus]
            ubV['BusVoltages',:,'Imag',bus] =  self.PowerFlowBounds['Vmax'][bus]
            lbV['BusVoltages',:,'Real',bus] = -self.PowerFlowBounds['Vmax'][bus]
            lbV['BusVoltages',:,'Imag',bus] = -self.PowerFlowBounds['Vmax'][bus]
          
        ubV["BusVoltages",:,"Imag",0] = 0.
        lbV["BusVoltages",:,"Imag",0] = 0.
        

        ######  EMBBED INITIAL CONDITIONS   #######
        if (self._hasStates == True):
            print "Initial Condition embedding"
            lbV['States',0] = x0 
            ubV['States',0] = x0

                
        ###### PERIODIC CONSTRAINTS (IF REQUIRED)   #######
        #if (Periodic == False):
        #    lbg['Periodic'] = -inf
        #    ubg['Periodic'] =  inf
        
        EP = self._EP()
        if not(ExtParameters == []):
            EP['ExtParameters'] = ExtParameters    
        EP['u0'] = u0
       
        self.OptDispatch.setInput(lbV,      "lbx")
        self.OptDispatch.setInput(ubV,      "ubx")
        self.OptDispatch.setInput(init,     "x0" )
        self.OptDispatch.setInput(lbg,      "lbg")
        self.OptDispatch.setInput(ubg,      "ubg")
        
        self.OptDispatch.setInput(EP,       "p")
        
        self.OptDispatch.solve()
        
        self.lbV = lbV
        self.ubV = ubV
        self.ubg = ubg
        self.lbg = lbg
        self.ep  = EP
        
        v_opt = self.VOptDispatch(self.OptDispatch.output("x"))
        
        success = int(self.OptDispatch.getStat('return_status') == 'Solve_Succeeded')
    
        return v_opt, success
    
    
    def Shift(self, Sol):
        
        SolShifted = self.VOptDispatch()
        
        for key in Sol.keys():
            Nelements = len(Sol[key])
            IndexTime =     [k for k in range( 1,Nelements ) ]
            IndexTimePlus = [k for k in range( Nelements-1 ) ]

            SolShifted[key,IndexTimePlus] = Sol[key,IndexTime]
            SolShifted[key,-1]            = Sol[key,-1]           
        return self.VOptDispatch(SolShifted)
    
    
    def Simulate(self, Sol, x0, u0):

        #To be replaced by a genuine simulation in the future...
        x0plus = self.x0(Sol['States',1])
        u0plus = self.u0(Sol['Inputs',0])
                
        return x0plus, u0plus
    
    
    
    def NMPCSimulation(self, x0 = [], u0 = [], ExtParameters = [], init = [], Simulation = 0):
        #####    NMPC Loop     #####
        NMPC = {'time': 0, 'success' : [], 'Traj' : []}
        Vstore = self.Vstore

        while (NMPC['time'] < Simulation):
            Sol, stats = self.DYNSolve(x0 = x0, u0 = u0, ExtParameters = ExtParameters, time = NMPC['time'], init = init)
            
            NMPC['success'].append(stats)
            NMPC['Traj'].append(Sol)
            
            Vstore[...,NMPC['time']] = Sol[...,0]
        
            init = self.Shift(Sol)    
            x0, u0 = self.Simulate(Sol,x0, u0)
            
            NMPC['time'] += 1
                               
        EP       = self._EP()                       
        EP['u0'] = u0
                                  
        self.LagrangeCost.setInput(Vstore,0)
        self.LagrangeCost.setInput(EP,1)
        self.LagrangeCost.evaluate()
        
        NMPC['LagrangeCost'] = self.LagrangeCost.output()
        return Vstore, NMPC    

    
    #Extract results    
    def ExtractInfo(self, v_opt, BusPower = True, PlantPower = True, TotalPower = True):
        
        self.SolutionInfo = {}
        
        
        Nstage = len(v_opt['Inputs'])
        
        NBus = self.NBus
        
        NLine     =  len(  self.Graph                )


        

        #DEFAULT EXTRACTION
        #Bus voltages (module and angles), Line Currents
        
        self.SolutionInfo['BusVoltagesModule'] = np.concatenate([np.array(np.sqrt(v_opt["BusVoltages",k,"Real",veccat]*v_opt["BusVoltages",k,"Real",veccat] + v_opt["BusVoltages",k,"Imag",veccat]*v_opt["BusVoltages",k,"Imag",veccat])).T for k in range(Nstage)],axis=0)
        self.SolutionInfo['BusVoltagesAngle'] = np.concatenate([np.array([180*math.atan2(v_opt["BusVoltages",k,"Imag",bus],v_opt["BusVoltages",k,"Real",bus])/pi for bus in range(NBus)]).reshape(NBus,1) for k in range(Nstage)], axis = 1).T
        
        LineCurrents_opt = []
        for k in range(Nstage):
            self.LineCurrents2Func.setInput(v_opt["BusVoltages",k])
            self.LineCurrents2Func.evaluate()
            LineCurrents_opt.append(sqrt(self.LineCurrents2Func.output()))
            
        
        self.SolutionInfo['LineCurrentsModule'] = np.concatenate(LineCurrents_opt,axis=1)
        
        #### Total Powers
        if (TotalPower == True):
            TotalPower = {}
            
            TotalPower['Load'] = 0
            for plant in [plant for plant in self.PlantList if (plant._Load == True)]:
                TotalPower['Load'] -= np.array(v_opt['Inputs',:,plant.label,'ActivePower'])
                   
            TotalPower['Injected'] = 0.
            for plant in [plant for plant in self.PlantList if not(plant._Load == True)]:     
                    if (plant.Directionality == 'Mono'):
                        TotalPower['Injected'] += np.array(v_opt['Inputs',:,plant.label,'Power'])
                    else:
                        TotalPower['Injected'] += np.array(v_opt['Inputs',:,plant.label,'Pdischarge']) - np.array(v_opt['Inputs',:,plant.label,'Pcharge'])        

            self.SolutionInfo['TotalPower'] = TotalPower             
        
        
        #Construct implicit values in the network 
        if (BusPower == True):
            
            BusActivePower = []
            BusReactivePower = []
            BusCurrentModule = []
            BusCurrentAngle = []
            for k in range(Nstage):
                self.BusActivePowerFunc.setInput(  v_opt['BusVoltages',k],0)
                self.BusReactivePowerFunc.setInput(v_opt['BusVoltages',k],0)
                self.BusCurrentsRealFunc.setInput( v_opt['BusVoltages',k],0)
                self.BusCurrentsImagFunc.setInput( v_opt['BusVoltages',k],0)
                self.BusActivePowerFunc.evaluate()
                self.BusReactivePowerFunc.evaluate()
                self.BusCurrentsRealFunc.evaluate()
                self.BusCurrentsImagFunc.evaluate()
            
                BusCurrentReal_k = np.array(self.BusCurrentsRealFunc.output())
                BusCurrentImag_k = np.array(self.BusCurrentsImagFunc.output())
                BusActivePower.append(       np.array(self.BusActivePowerFunc.output()).T)
                BusReactivePower.append(   np.array(self.BusReactivePowerFunc.output()).T)
                BusCurrentModule.append(   sqrt(BusCurrentReal_k**2 + BusCurrentImag_k**2).T )
                BusCurrentAngle.append( np.array([180*math.atan2(BusCurrentImag_k[bus],BusCurrentReal_k[bus])/pi for bus in range(NBus)  ]).reshape(1,6) )
                
            self.SolutionInfo['BusActivePower']   = np.concatenate(  BusActivePower,   axis=0)
            self.SolutionInfo['BusReactivePower'] = np.concatenate(  BusReactivePower, axis=0)
            self.SolutionInfo['BusCurrentModule'] = np.concatenate(  BusCurrentModule, axis=0)
            self.SolutionInfo['BusCurrentAngle']  = np.concatenate(  BusCurrentAngle,  axis=0)
        
        if (PlantPower == True):
            PlantActivePowerDictionary = {}
            PlantReactivePowerDictionary = {}
            CosPhiDictionary = {}

            for plant in self.PlantList:
                #Data = self.Plants[key]
                #Nplant = np.size(Data,axis = 0)
                PlantActivePower = []
                PlantReactivePower = []
                CosPhi = []
                TanPhi = []

                Bus = plant.Bus
                CurrentReal =    veccat(v_opt['Inputs',:,plant.label,'CurrentReal'])
                CurrentImag =    veccat(v_opt['Inputs',:,plant.label,'CurrentImag'])
                BusVoltageReal = veccat(v_opt['BusVoltages',:,'Real',Bus])
                BusVoltageImag = veccat(v_opt['BusVoltages',:,'Imag',Bus])
                
                PlantActivePower_plant = np.array(BusVoltageReal*CurrentReal + BusVoltageImag*CurrentImag)
                PlantReactivePower_plant = np.array(BusVoltageImag*CurrentReal - BusVoltageReal*CurrentImag)
                PlantApparentPower_plant = sqrt(PlantReactivePower_plant**2 + PlantActivePower_plant**2)
                
                PlantActivePower.append(PlantActivePower_plant)
                PlantReactivePower.append(PlantReactivePower_plant)
                CosPhi.append(PlantActivePower_plant/PlantApparentPower_plant)
                    
                    
                PlantActivePowerDictionary[plant.label]   =  PlantActivePower
                PlantReactivePowerDictionary[plant.label] =  PlantReactivePower
                CosPhiDictionary[plant.label]             =  CosPhi
                
            self.SolutionInfo['PlantActivePower']         =  PlantActivePowerDictionary
            self.SolutionInfo['PlantReactivePower']       =  PlantReactivePowerDictionary
            self.SolutionInfo['PlantCosPhi']              =  CosPhiDictionary



    
    #### RESULT PLOTTING #####
    

    def DYNSolvePlot(self, v_opt, NMPC = False, dt = 1, Path = [], LW = 1, Show = True):
        
        SavedFigs = []
        
        def SaveFig(Path,Name):
            if not(Path == []):
                SavedFigs.append(Name)
                plt.tight_layout()
                plt.savefig(Path+'/'+Name+'.eps',format='eps', facecolor='w', edgecolor='k',bbox_inches='tight')
                #plt.close()
            return SavedFigs
        
        
        if (len(v_opt['States']) < 2):
            print "Plotting warning: no time sequence available, run with .Nsim > 1. Plotting not proceeding."
            return
        
        #Nstage = self.TimeSetup['Horizon']

        NBus = self.NBus
        NLine     =  len(    self.Graph   )


        
        #####  Prepares time grids #####

        #dt = self.TimeSetup['dt']
        time = {}
        for key in ['States','Inputs']:
            time[key] = np.array([k*dt for k in range(len(v_opt[key]))]).T
        
        # construct a list of the plants (excluding the loads)
        PlantList = [plant for plant in self.PlantList if not(plant._Load == True)]
        SizeSubplt    = np.ceil(sqrt(len(PlantList)))
        SizeSubpltAll = np.ceil(sqrt(len(self.PlantList)))
        
        ##### plot plot plot #####
        
        plt.figure(1)
        plt.subplot(2,3,1)
        plt.hold('on')
        for k in range(NBus):
            plt.step(time['Inputs'],self.SolutionInfo['BusVoltagesModule'][:,k],where = 'post', label = str(k), linewidth = LW)
        plt.ylabel('kV')
        #plt.xlabel('time (s)')
        plt.title("Voltages, |.|")
        plt.grid()
        
        plt.subplot(2,3,2)
        plt.hold('on')
        for k in range(NBus):
            plt.step(time['Inputs'],self.SolutionInfo['BusVoltagesAngle'][:,k],where = 'post', label = str(k), linewidth = LW)
        plt.ylabel('deg')
        #plt.xlabel('time (s)')
        plt.title("Voltage, angle")
        plt.grid()
        
        plt.subplot(2,3,5)
        plt.hold('on')
        for k in range(NBus):
            plt.step(time['Inputs'],1e-3*self.SolutionInfo['BusActivePower'][:,k],where = 'post', label = str(k), linewidth = LW)
        plt.ylabel('GW')
        plt.xlabel('time (s)')
        plt.title("Active power")
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        plt.grid()
        
        plt.subplot(2,3,4)
        plt.hold('on')
        for k in range(NBus):
            plt.step(time['Inputs'],1e-3*self.SolutionInfo['BusReactivePower'][:,k],where = 'post', label = str(k), linewidth = LW)
        plt.ylabel('GW')
        plt.xlabel('time (s)')
        plt.title("Reactive power")
        plt.grid()
        
        plt.subplot(2,3,3)
        plt.hold('on')
        for k in range(NBus):
            plt.step(time['Inputs'],self.SolutionInfo['BusCurrentModule'][:,k],where = 'post', label = str(k), linewidth = LW)
        plt.ylabel('kA')
        plt.xlabel('time (s)')
        plt.title('Current, |.|')
        #plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        plt.grid()
        
        SaveFig(Path,'Grid')
        
        
                
        plt.figure(2)
        for k in range(NLine):
            plt.step(time['Inputs'],self.SolutionInfo['LineCurrentsModule'][k,:],where = 'post', label = str(self.Graph[k][0])+'-'+str(self.Graph[k][1]), linewidth = LW)
        
        plt.xlabel('time (s)')
        plt.ylabel('kA')
        plt.title("Lines current |.|")
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        plt.grid()

        SaveFig(Path,'Lines')

                
        plt.figure(7)
        fig = 1
        for plant in self.PlantList:            
            plt.subplot(SizeSubpltAll,SizeSubpltAll,fig)
            plt.step(time['Inputs'],1e-3*self.SolutionInfo['PlantActivePower'][plant.label][0], color = 'k', label = 'Act. Power', linewidth = LW)
            plt.step(time['Inputs'],1e-3*self.SolutionInfo['PlantReactivePower'][plant.label][0], color = 'r', label = 'React. Power', linewidth = LW)
            
            fig += 1
            plt.ylabel('GW')
            plt.xlabel('time (s)')
            plt.title(str(plant.label)+', bus '+str(plant.Bus))
            plt.grid()
            
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

        SaveFig(Path,'BusPower')

        
        plt.figure(8)
        plt.hold('on')
        for plant in self.PlantList:
            Power = []
            if   plant._Load == True:
                Power = -np.array(v_opt['Inputs',:,plant.label,'ActivePower'])
            elif plant.Directionality == 'Bi':
                Power = np.array(v_opt['Inputs',:,plant.label,'Pdischarge']) - np.array(v_opt['Inputs',:,plant.label,'Pcharge'])
            elif plant.Directionality == 'Mono':
                Power = v_opt['Inputs',:,plant.label,'Power']
            else:
                print "Warning: plant power unidentified, not plotting"
            
            if len(Power)>0:
                plt.step(time['Inputs'],1e-3*np.array(Power), label = plant.label,where = 'post', linewidth = LW)
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        plt.xlabel('time (s)')
        plt.ylabel('GW')
        plt.title('Plant power')
        plt.grid()
        
        SaveFig(Path,'PlantPower')
 
        

        plt.figure(9)
        plt.step(time['Inputs'],100*self.SolutionInfo['TotalPower']['Load']/self.SolutionInfo['TotalPower']['Injected'], label = plant.label,where = 'post', color = 'k', linewidth = LW)
        plt.xlabel('time (s)')
        plt.ylabel('%')
        plt.title('Transmission efficiency')
        plt.grid()
        
        SaveFig(Path,'GridEfficiency')

        
        plt.figure(10)
        fig = 0
        Nsubplot = 0
        UnitDic = {'h': 'm', 'W_error': 'm/s','E': 'MJ'}
        for plant in PlantList:
            if hasattr(plant,'_Shoot'):
                Nsubplot += len(plant.States.keys())
        Nsubplot = ceil(sqrt(Nsubplot))
        for plant in PlantList:
            if hasattr(plant,'_Shoot'):
                for key in plant.States.keys():
                    plt.subplot(Nsubplot,Nsubplot,fig)
                    plt.step(time['States'],np.array(v_opt['States',:,plant.label,key]),color = 'k', linewidth = LW)
                    plt.title(plant.label+', state: '+key)
                    plt.xlabel('time (s)')
                    plt.ylabel(UnitDic[key])
                    plt.grid()
                    fig += 1
                    
        SaveFig(Path,'PlantStates')


        #Plant Detail
        fig = 11
        for plant in PlantList:
            if hasattr(plant,'_additionalInputs') or hasattr(plant,'_Shoot'):
                plt.figure(fig)
                plt.title(plant.label)
                subPltNum = 2
                if hasattr(plant,'_additionalInputs'):
                    subPltNum += len(plant._additionalInputs)
                if hasattr(plant,'_Shoot'):
                    subPltNum += len(plant.States.keys())
                subPltNum = np.ceil(sqrt(subPltNum))
            
                #Plot current
                plt.subplot(subPltNum,subPltNum,0)
                plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,'CurrentReal']),where = 'post',label = 'Real Current')
                plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,'CurrentImag']),where = 'post',label = 'Complex Current')
                plt.title('Current')
                plt.legend()
                plt.grid()
                
                plt.subplot(subPltNum,subPltNum,1)
                plt.title('Power')
                if plant.Directionality == 'Mono':
                    plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,'Power']),where = 'post',label = 'Power') 
                else:
                    plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,   'Pcharge']),where = 'post',label = 'Pcharge')
                    plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,'Pdischarge']),where = 'post',label = 'Pdischarge')
                plt.legend()    
                plt.grid()
                
                subplt = 2                
                if hasattr(plant,'_additionalInputs'):
                    for key in plant._additionalInputs:
                        subplt += 1
                        plt.subplot(subPltNum,subPltNum,subplt)
                        plt.step(time['Inputs'],np.array(v_opt['Inputs',:,plant.label,key]),where = 'post',label = key) 
                        plt.legend()
                        plt.grid()
                if hasattr(plant,'_Shoot'):
                    for key in plant.States.keys():
                        subplt += 1
                        plt.subplot(subPltNum,subPltNum,subplt)
                
                        plt.plot(time['States'],np.array(v_opt['States',:,plant.label,key]),label = key) 
                        plt.legend()#bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
                        plt.grid()
                        
                fig += 1
                
        if (Show == True):        
            plt.show()
        else:
            plt.close('all')
            
        return SavedFigs
    
    def DYNSolvePlotCompare(self, v_opt, NMPC = False, label = '', marker = '', linewidth = 1):
    
        Nstage = self.TimeSetup['Horizon']

        NBus = self.NBus
        NLine     =  len(    self.Graph   )

        NLoad     =  len(  self.Plants['Load']       )
        NStorage  =  len(  self.Plants['Storage']    )
        NWind     =  len(  self.Plants['Wind']       )
        NThermal  =  len(  self.Plants['Thermal']    )
        NHydro    =  len(  self.Plants['Hydro']      )
        
        

        dt = self.TimeSetup['dt']
        time = {}
        for key in ['States','Inputs']:
            time[key] = np.array([k*dt for k in range(len(v_opt[key]))]).T
        

        plt.figure(3)
        fig = 1
        for key in ['Wind','Hydro','Thermal','Storage']:
            plt.subplot(2,2,fig)
            plt.hold('on')
            plt.title(key+' Power')
            for plant in range(eval('N'+key)):
                if (self.Plants[key][plant]['Directionality'] == 'mono'):
                    plt.step(time['Inputs'],np.array(v_opt['Inputs',:,key,plant,'Power']),where = 'post',label = label, linewidth = linewidth)
                    if (key == 'Wind'):
                        CPmax = self.Plants['Wind'][plant]['CPmax']
                        A     = self.Plants['Wind'][plant]['A']
                        PWind = 0.5*rho_air*A*CPmax*np.array(v_opt['States',:,'Wind',plant,'WindSpeed'])**3                
                        plt.step(time['States'],PWind,color='k',where = 'post')
                        plt.hlines(self.Plants['Wind'][plant]['Pmax'], time['States'][0], time['States'][-1], colors='r', linestyles='dashed')
                else:
                    plt.step(time['Inputs'],np.array(v_opt['Inputs',:,key,plant,'Pdischarge'])-np.array(v_opt['Inputs',:,key,plant,'Pcharge']),where = 'post',label = label, linewidth = linewidth) 
        
            fig += 1
            
        #########    

        plt.figure(4)
        plt.subplot(2,2,1)
        plt.hold('on')
        for plant in range(NHydro):
            plt.plot(time['States'],np.array(v_opt['States',:,'Hydro',plant,'WaterHeight']),label = label, marker = marker, linewidth = linewidth)
        plt.title('Water Height')
        
        plt.subplot(2,2,2)
        plt.hold('on')
        for plant in range(NStorage):
            plt.plot(time['States'],np.array(v_opt['States',:,'Storage',plant,'Energy']),label = label, marker = marker, linewidth = linewidth)
        plt.title('Stored Energy')
        
        plt.subplot(2,2,3)
        plt.hold('on')
        for plant in range(NHydro):
            plt.step(time['Inputs'],sqrt(np.array(v_opt['Inputs',:,'Hydro',plant,'Pcharge'])*np.array(v_opt['Inputs',:,'Hydro',plant,'Pdischarge'])),where = 'post',label = label, marker = marker, linewidth = linewidth)
        plt.title('Hydro Complementarity')
        
        plt.subplot(2,2,4)
        plt.hold('on')
        for plant in range(NStorage):
            plt.step(time['Inputs'],sqrt(np.array(v_opt['Inputs',:,'Storage',plant,'Pcharge'])*np.array(v_opt['Inputs',:,'Storage',plant,'Pdischarge'])),where = 'post',label = label, marker = marker, linewidth = linewidth)
        plt.title('Storage Complementarity')
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        
        plt.figure(8)
        for plant in range(NWind):
            plt.step(time['States'],veccat(v_opt['States',:,'Wind',plant,'WindSpeed']),label=label, marker = marker, linewidth = linewidth)
        plt.title('Wind Speed m/s')
                

 
