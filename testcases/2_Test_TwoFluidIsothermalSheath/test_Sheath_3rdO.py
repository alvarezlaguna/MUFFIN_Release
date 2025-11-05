#!/usr/bin/env python3
import sys
from mpi4py import MPI
sys.path.append(sys.path[0] + '/..')
import muffin				## module in C++ doing the simulation
import glob				## module for unix pathname pattern expansioni
import matplotlib			## Add these two lines to show results
matplotlib.use('TkAgg') 		## It uses tk to not freeze
import matplotlib.pyplot as plt 	## plots 
import os 				## operating system interface
import numpy as np			## Mathematic operations
import shutil				## removing the old result directory

comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
world_size = comm.size 

# Dimensionless parameters for the simulation
MassRatio        = 1e5 		# mi/me
TemperatureRatio = 0.025 	# Ti/Te
Debye 		     = 1.		# L/lambda (I usually take 1 and take it as reference length. Then you need to adjust the size of the domain)
CollIons         = 0.		# nu_in/t_0 (The reference time is taken as t_0 = L_0/u_B. Note that if the Debye = 1, the L_0 = Debye)
CollElectrons    = 0.           # nu_en/t_0
### Note that the simulation computes the ionization frequency according to the flux on the wall, so not need to change 'ionizConst' does not need to be changed 

## Set the options
nbCells = 10000
InitialField = [0 for i in range(4*nbCells)]
# Set Density of Electrons
for iCell in range(nbCells):
       InitialField[0*nbCells + iCell] = 1.
# Set Momentum of Electrons
for iCell in range(nbCells):
       InitialField[1*nbCells + iCell] = 0.
# Set Density of Ions
for iCell in range(nbCells):
       InitialField[2*nbCells + iCell] = 1.
# Set Momentum of Ions
for iCell in range(nbCells):
       InitialField[3*nbCells + iCell] = 0.

# Set Mesh
length = 100.		# Length of the domain in Debye lengths
mesh = np.linspace(0., length, nbCells+1)

# Set Phi Initial
PhiInitial = [0 for i in range(nbCells)]

options = {
	'nbEqs':4,
	'nbFluids':2,
	'nbCells':[nbCells],
	'geometry':"1D",
	'mesh':mesh,
 	'length':length,
	#'stopCondition':{'type':"nbSteps",'value':100},
	'stopCondition':{'type':"Residual",'value':-7},
	'Inlet':{'type':"TwoFluidIsothermalSheathHagelaar",},
	'Outlet':{'type':"TwoFluidIsothermalSheathHagelaar",},
	'PhysicalModel':{
		'type':"MultiFluidIsothermal1D",
		'soundSpeed': [np.sqrt(MassRatio), np.sqrt(TemperatureRatio)],
		},
	'SourceTerm':{
		'type':"TwoFluidIsothermal1D",
		'massRatio':MassRatio,
		'DebyeLength':Debye,
		'ionizConst':1.,
		'epsIoniz':8.72,
		'CollIons':CollIons,
		'CollElectrons':CollElectrons,
		'PhiIn':0.,
		'PhiOut':0.,
		'PhiInitial':PhiInitial,
		},
	'CFL':0.9,
	'saveRate':1000,
	'resultDir':"./test_Sheath",

	'limiter':"thirdOrder", #"Venkatakrishnan", #"ospre", #"thirdOrder",
	'reconstruction':"TVD2ndOrder1D",#"1stOrder"
	'fluxScheme':"RoeMultiFluidIsothermal1D",#"LaxFriedrich",
    'timeScheme':"TVDRK3",#"ForwardEuler",#"TVDRK3",
	'initialField':InitialField,
}

## Run simulation
if os.path.exists(options['resultDir']):
    shutil.rmtree(options['resultDir']) # Removes the olde results
comm.barrier()
muffin.Solver(options)
sys.stdout.flush()


