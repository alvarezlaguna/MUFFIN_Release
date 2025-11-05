import sys
from mpi4py import MPI
sys.path.append(sys.path[0] + '/../../PythonTools/')
import plotResults
sys.path.append('../../build')                     ## We add this to include the module muffin
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

# Set the density
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
        #'stopCondition':{'type':"Residual",'value':-12},
        'stopCondition':{'type':"nbSteps",'value':200},
        'saveRate':5,
        'resultDir':"./Results_SodTube_FirstOrder_Roe",
        'limiter':"ospre",#"Venkatakrishnan",#"vanAlbada",#"ospre",#"thirdOrder"
        'reconstruction':"TVD2ndOrder1D",#"1stOrder","TVD2ndOrder1D"
        'fluxScheme':"RoeEuler1D",#"RoeEuler1D", "LaxFriedrich","HLLEuler1D"
        'timeScheme':"ForwardEuler",#"TVDRK3",
        'initialField':InitialField,
}

## Remove previous simulation
if os.path.exists(options['resultDir']):
    shutil.rmtree(options['resultDir']) # Removes the old results

## Run simulation
#muffin.Solver(options)
simu=muffin.Simulation1D("bonjour")
simu.setup(options)
a=simu.simulate()
print(a)

## Copy input file inside the Result folder
os.system("cp "+os.path.basename(__file__)+" "+options['resultDir'])

