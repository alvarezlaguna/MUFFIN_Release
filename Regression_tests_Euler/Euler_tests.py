import sys
from mpi4py import MPI
import muffin                                   ## module in C++ doing the simulation
import glob                                     ## module for unix pathname pattern expansioni
import matplotlib.pyplot as plt                 ## plots
import os                                       ## operating system interface
import numpy as np                              ## Mathematic operations:w
import shutil                                   ## Remove folder


## Set the options
nbCells_half = 100
nbCells = int(2*nbCells_half)
nbEqs   = 3
InitialField = [0 for i in range(3*nbCells)]
length  = 1.

# Set Mesh
mesh = [0 for i in range(nbCells + 1)]
for iCell in range(nbCells + 1):
       delta_X = length/nbCells
       mesh[iCell] = iCell*delta_X

# # Set the density
for iCell in range(nbCells_half):
        InitialField[0*nbCells + iCell] = 1
for iCell in range(nbCells_half, nbCells):
        InitialField[0*nbCells + iCell] = 0.125
# Set the momentum
for iCell in range(nbCells):
        InitialField[1*nbCells + iCell] = 0
# Set the total energy
for iCell in range(nbCells_half):
        InitialField[2*nbCells + iCell] = 1./(5./3. - 1.)
for iCell in range(nbCells_half, nbCells):
        InitialField[2*nbCells + iCell] = 0.1/(5./3. - 1.)



def Euler_first_order(scheme):

        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                'type':"EulerEq1D",
                'gamma':5./3.,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./1_FirstOrder_Roe",
                'limiter':"vanAlbada",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"1stOrder",#"TVD2ndOrder1D"
                'fluxScheme':scheme,#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"ForwardEuler",#"TVDRK3",
                'initialField':InitialField,
        }

        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]

def Euler_TVDRK3():

        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                'type':"EulerEq1D",
                'gamma':5./3.,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./2_FirstOrder_Roe_TVDRK3",
                'limiter':"vanAlbada",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"1stOrder",#"TVD2ndOrder1D"
                'fluxScheme':"RoeEuler1D",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"TVDRK3",#"TVDRK3",
                'initialField':InitialField,
        }

        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]

def Euler_SSPRK3():

        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                'type':"EulerEq1D",
                'gamma':5./3.,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./3_FirstOrder_Roe_SSPRK3",
                'limiter':"vanAlbada",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"1stOrder",#"TVD2ndOrder1D"
                'fluxScheme':"RoeEuler1D",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"SSPRK3",#"TVDRK3",
                'initialField':InitialField,
        }

        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]

def Euler_SSPRK5():

        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                'type':"EulerEq1D",
                'gamma':5./3.,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./3_FirstOrder_Roe_SSPRK5",
                'limiter':"vanAlbada",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"1stOrder",#"TVD2ndOrder1D"
                'fluxScheme':"RoeEuler1D",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"SSPRK5",#"TVDRK3",
                'initialField':InitialField,
        }

        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]

def Euler_Roe_second_order(limiter):

        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                'type':"EulerEq1D",
                'gamma':5./3.,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./Results_SodTube_FirstOrder_Roe",
                'limiter':limiter,#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"TVD2ndOrder1D",
                'fluxScheme':"RoeEuler1D",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"ForwardEuler",#"TVDRK3",
                'initialField':InitialField,
        }


        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]

def Euler_first_order_pythonModel():


        def flux_function(u, m_physFlux):
                
                gamma  = 5./3.
                rho    = u[0]
                m      = u[1]
                e      = u[2]
                
                m_physFlux[0] = m
                m_physFlux[1] = m*m/rho*(3 - gamma)/2 + (gamma - 1)*e
                m_physFlux[2] = e*m/rho*gamma - (gamma - 1)/2*m**3/rho**2


        def maxEigen_function(u):
                gamma  = 5./3.
                rho    = u[0]
                m      = u[1]
                e      = u[2]
                v      = m/rho
                soundSpeed = np.sqrt(gamma*(gamma - 1)*(e - m*m/(2*rho))/rho)
    
                m_eigenvals_1 = v - soundSpeed
                m_eigenvals_2 = v
                m_eigenvals_3 = v + soundSpeed

                return np.max([np.abs(m_eigenvals_1), np.abs(m_eigenvals_2), np.abs(m_eigenvals_3)])
                


        options = {
                'nbEqs':nbEqs,
                'nbFluids':1,
                'nbCells':[nbCells],
                'geometry':"1D",
                'mesh':mesh,
                'length':length,
                'Inlet':{'type':"Neumann",},
                'Outlet':{'type':"Neumann",},
                'PhysicalModel':{
                        'type':"PythonPhysicalModel",
                        'flux_function':flux_function,
                        'maxEigen_function': maxEigen_function,
                },
                'SourceTerm':{'type':"NullSourceTerm",},
                'CFL':0.3,
                'stopCondition':{'type':"nbSteps",'value':100},
                'saveRate':100,
                'resultDir':"./1_FirstOrder_Roe",
                'limiter':"vanAlbada",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
                'reconstruction':"1stOrder",#"TVD2ndOrder1D"
                'fluxScheme':"LaxFriedrich",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
                'timeScheme':"ForwardEuler",#"TVDRK3",
                'initialField':InitialField,
        }

        ## Run simulation
        Residual = muffin.Solver(options)

        ## Copy input file inside the Result folder
        os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

        return Residual, options["resultDir"]









