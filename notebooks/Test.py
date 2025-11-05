import sys
from mpi4py import MPI
sys.path.append(sys.path[0] + '/../PythonTools/')
import plotResults 
sys.path.append(sys.path[0] + '/../build/')     ## We add this to include the module muffin
print(sys.path[0] + '/../build/')
import muffin                                   ## module in C++ doing the simulation
import glob                                     ## module for unix pathname pattern expansioni
import matplotlib.pyplot as plt                 ## plots
import os                                       ## operating system interface
import numpy as np                              ## Mathematic operations

# Set Mesh

length  = 1.
nbCells = 100 

mesh = [0 for i in range(nbCells + 1)]
for iCell in range(nbCells + 1):
    delta_X = length/nbCells 
    mesh[iCell] = iCell*delta_X

nbEqs   = 1
InitialField = [0 for i in range(nbEqs*nbCells)]

options = {
    'nbEqs':nbEqs,
    'nbFluids':1,
    'nbCells':[nbCells],
    'geometry':"1D",
    'mesh':mesh,
    'length':length,
    'stopCondition':{'type':"nbSteps",'value':2000},
    #'stopCondition':{'type':"Residual",'value':-12},
    'Inlet':{'type':"Dirichlet", 'value':[1.],},
    'Outlet':{'type':"Neumann",},
    'PhysicalModel':{
        'type':"AdvectionEq1D",
        'A':1.,
    },
    'SourceTerm':{'type':"NullSourceTerm",}, 
    'CFL':0.1,
    'saveRate':20,
    'resultDir':"./Results_advection",
    'limiter':"thirdOrder",#"thirdOrder",
    'reconstruction':"1stOrder",#"TVD2ndOrder1D",
    'fluxScheme':"LaxFriedrich",
    'timeScheme':"ForwardEuler",#"TVDRK3",
    'initialField':InitialField,
}

## Run simulation
if os.path.exists(options['resultDir']):
    shutil.rmtree(options['resultDir']) # Removes the old results

muffin.Solver(options)

sys.stdout.flush()

