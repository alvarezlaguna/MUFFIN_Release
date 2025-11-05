//
//  Simulate1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::min_element, std::max_element
#include <numeric>

#include "Simulation1D.hpp"
#include "MeshData/MeshData.hpp"
#include "MeshData/Cell1D.hpp"
#include "Options.hpp"

#include "FluxScheme/FluxScheme.hpp"

#include "SpaceReconstructor/SpaceReconstructor.hpp"
#include "SpaceReconstructor/ReconstructorFactory.hpp"
#include "SpaceReconstructor/Limiter/Limiter.hpp"
#include "SpaceReconstructor/Limiter/LimiterFactory.hpp"

#include "FluxScheme/FluxSchemeFactory.hpp"
#include "SpaceMethod/FVM1D.hpp"
#include "TimeMethod/TimeMethodFactory.hpp"
#include "PhysicalModel/PhysicalModelFactory.hpp"
#include "DataWriter/DataWriterFactory.hpp"

#include "CollisionalData/CollisionalData.hpp"
#include "CollisionalData/Mixture.hpp"
#include "CollisionalData/MixtureFactory.hpp"
#include "CollisionalData/SingleIonMixture.hpp"

#include "DiffusiveFlux/DiffusiveFlux.hpp"
#include "DiffusiveFlux/DiffusiveFluxFactory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <mpi.h>

#include "SourceTerm/InputFilesProvider.hpp"

namespace py = pybind11;
using namespace std;
using namespace Parameters;

// -------------------------------------------------------------------------- //

/// function that generates the mesh data

void Simulation1D::setup() {
    
    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::locale::global(std::locale::classic());

    if(rank == MPI_WRITER)
        py::print("Setting-up Simulation class:\t ",this->getName(),"\n");

    // Generate the general mesh data from mother class
    Simulation::setup();
    
    // Generating extra data for 1D
    MeshData& md = MeshData::getInstance();
    // Create the data
    md.createData<Cell1D>("Cells", NBCELLS);
    md.create2DData<double>("Cells_py_CC", NBEQS, NBCELLS+2); 
    md.create2DData<double>("Cells_py_L", NBEQS,NBCELLS);
    md.create2DData<double>("Cells_py_R", NBEQS,NBCELLS);
    md.createData<double>("x", NBCELLS);                // x-Coordinate of the mesh
    md.create2DData<double>("source",NBEQS, NBCELLS);     // Source term
    md.createData<CellDataRef>("boundaries",2);        // 2 boundary conditions needed: Inlet and Outlet
    
    // Create the data
    vector<Cell1D>& cells = md.getData<Cell1D>("Cells");
    py::array_t<double>& py_cells_CC = md.get2DData<double>("Cells_py_CC");
    py::array_t<double>& py_cells_L = md.get2DData<double>("Cells_py_L");
    py::array_t<double>& py_cells_R = md.get2DData<double>("Cells_py_R");
    for (int iCell = 0; iCell < Parameters::NBCELLS; iCell++){
        cells[iCell]=Cell1D(CellDataRef(&py_cells_CC,iCell+1), CellDataRef(&py_cells_L,iCell), CellDataRef(&py_cells_R,iCell));
    }
    vector<CellDataRef>& boundaries = md.getData<CellDataRef>("boundaries");
    boundaries[0]=CellDataRef(&py_cells_CC,0);
    boundaries[1]=CellDataRef(&py_cells_CC,NBCELLS+1);
    // Set-up the intial field
    setInitialField();

    // Initialize the 1D mesh
    vector<double>& x = MeshData::getInstance().getData<double>("x");

    const int n1 = NBCELLS - 1;
    // First cell
    double delta_X  = MESH[1] - MESH[0];
    x[0]            = (MESH[1] + MESH[0])/2.;
    cells[0].cellID = 0;
    cells[0].dx     = delta_X;
    
    // Inner cells
    for (int i = 1; i < n1; i++){   // Set-up the mesh
        delta_X = (MESH[i + 1] - MESH[i]);
        
        // Set the coordinate
        x[i] = (MESH[i + 1] + MESH[i])/2;
        
        // Set the ID of the cell
        cells[i].cellID = i;
        
        // Set the dx of the cell
        cells[i].dx     = delta_X;
    }
    // Last cell
    delta_X          = MESH[NBCELLS] - MESH[n1];
    x[n1]            = (MESH[NBCELLS] + MESH[n1])/2;
    cells[n1].cellID = n1;
    cells[n1].dx     = delta_X;
    
    MESH.clear(); // deallocate the mesh vector
    
    // Setup the mixture
    auto mixType   = CollisionalData::getInstance().getData<string>("mixtureType");
    m_mixture      = MixtureFactory::CreateMixture(mixType);
    m_mixture->setup();
    
    auto diffFlux   = Options::getInstance().getData<string>("DiffusionType");
    m_diffFlux      = DiffusiveFluxFactory::CreateDiffusiveFlux(diffFlux);
    m_diffFlux->setMixture(m_mixture.get());
    m_diffFlux->setup();
    
    // Setup the physical model
    m_physicalModel = PhysicalModelFactory::CreatePhysicalModel(PHYSICALMODELNAME);
    m_physicalModel->setMixture(m_mixture.get());
    m_physicalModel->setDiffusiveFlux(m_diffFlux.get());
    
    // Setup the space method
    m_spaceMethod.reset(new FVM1D("FVM1D"));
    m_spaceMethod->setPhysicalModel(m_physicalModel.get());
    m_spaceMethod->setup();

    // Setup Time method
    m_timeMethod = TimeMethodFactory::CreateTimeMethod(TIMESCHEMENAME);
    m_timeMethod->setSpaceMethod(m_spaceMethod.get());
    m_timeMethod->setDiffusiveFlux(m_diffFlux.get());
    m_timeMethod->setup();
    
    // Setup Data Writer
    m_dataWriter = DataWriterFactory::CreateDataWriter(DATAWRITER);
    
    // Create result directory
    if (rank == 0){
        py::object os        = py::module::import("os");
        py::object makedirs  = os.attr("makedirs");
        py::object path      = os.attr("path");
        py::object pathExits = path.attr("exists");
        py::bool_ check      = pathExits(RESULTDIRECTORY);
        if (!check) {
            makedirs(RESULTDIRECTORY);
        }
    }
    if(rank == MPI_WRITER)
        py::print("Writing results in :\t",RESULTDIRECTORY,"\n");

}


// -------------------------------------------------------------------------- //

/// function that cleans up the mesh data

void Simulation1D::unsetup() {

    // unsetup the mother class
    Simulation::unsetup();
    MeshData& md = MeshData::getInstance();
    md.deleteData<Cell1D>("Cells");     // Data stored in the cells
    md.deleteData<double>("x");         //  x-Coordinate of the mesh
    md.delete2DData<double>("source");    // Source term
    md.deleteData<double>("boundaries");// boundary conditions
    
    // release the pointer
    m_spaceMethod.release();
    m_timeMethod.release(); 
}

// -------------------------------------------------------------------------- //

/// function that sets the initial field

void Simulation1D::setInitialField()
{
    //const int sizeU = NBCELLS*NBEQS;     // size of solution array
    // Initialize the 1D solution to a constant value
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    for(int iEq = 0; iEq < NBEQS; iEq++){
        for (int iCell = 0; iCell < NBCELLS; iCell++){
            cells[iCell].uCC[iEq] = INITIALFIELD[iEq*NBCELLS + iCell];
            cells[iCell].uL[iEq] = cells[iCell].uCC[iEq]; //We initialize the left and right states to the initial field
            cells[iCell].uR[iEq] = cells[iCell].uCC[iEq];

        }
    }
    
    INITIALFIELD.clear(); // deallocate the initial field vector
}

// -------------------------------------------------------------------------- //

/// function that writes data
void Simulation1D::writeData(const int iter)
{
    m_dataWriter->writeData(iter);

}

void Simulation1D::computeNorm(vector<double>& norm, vector<double>& normGlobal)
{
    
    // MPI excahnge of data
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // Compute the norm and send to all the processors
    MPI_Allreduce(norm.data(), normGlobal.data(), norm.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // Compute Norm
    for(int iEq = 0; iEq < NBEQS; ++iEq){
        if(normGlobal[iEq] != 0 ){ //check that the residual is not -inf
            normGlobal[iEq] = log10(sqrt(normGlobal[iEq]));
        }
        else{ normGlobal[iEq] = 0.; } //when the residual is 0, we set the log to 0 as well
    }
}

// -------------------------------------------------------------------------- //

/// Main loop of the 1D Finite Volume Simulation

py::array_t<double> Simulation1D::simulate(){
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    physTime = INITIAL_TIME;
    // Add a flag to see if the final time stop condition is true
    bool final_time_cond = (FINAL_TIME == 0) ? false : true;
    
    
    unsigned int iter = m_timeMethod->getIter();
    vector<double>& norm = m_timeMethod->getNorm();
    double meanNorm = EPS + 1; // We initialize to a value that is larger than the norm
    
    vector<double> norm_global(NBEQS);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // needed to store the norm and time in a buffer to use py::print
    char buffer_norm[100];
    char buffer_time[100];
    
    if(final_time_cond == false){
        while ((meanNorm > EPS) && iter < NBSTEPS) { // iterate till convergence

            // Time discretization
            m_timeMethod->takeStep();
            
            iter = m_timeMethod->getIter();
            norm = m_timeMethod->getNorm();
            computeNorm(norm, norm_global);
            meanNorm = 0;

            if (iter % SAVERATE == 0 || iter == 1){  writeData(iter);} //write data

            sprintf(buffer_norm, "%0.5lf",norm_global[0]);
            meanNorm += norm_global[0];

            for(int iEq =1; iEq<NBEQS; iEq++){
                sprintf(buffer_norm + strlen(buffer_norm), "\t%0.5lf",norm_global[iEq]);
                meanNorm += norm_global[iEq];
            }
            sprintf(buffer_time, "%0.3e",physTime);
            
            if (my_rank == MPI_WRITER){
                if (iter % PRINTRATE == 0){
                    py::print( "Iter[",iter,"] \t Res = [", buffer_norm , "]\t PhysTime = ", buffer_time);
                }
            }
            // Compute the mean of the norm to stop simulation
            meanNorm = meanNorm/NBEQS;        
        }
    }
    else{
        while (physTime < FINAL_TIME) { // iterate till convergence

            // Time discretization
            m_timeMethod->takeStep();
            
            iter = m_timeMethod->getIter();
            norm = m_timeMethod->getNorm();
            computeNorm(norm, norm_global);
            meanNorm = 0;

            if ( (SAVERATE!=0 && iter % SAVERATE == 0) || iter == 1){  writeData(iter);} //write data

            sprintf(buffer_norm, "%0.5lf",norm_global[0]);
            meanNorm += norm_global[0];

            for(int iEq =1; iEq<NBEQS; iEq++){
                sprintf(buffer_norm + strlen(buffer_norm), "\t%0.5lf",norm_global[iEq]);
                meanNorm += norm_global[iEq];
            }
            sprintf(buffer_time, "%0.3e",physTime);
            
            if (my_rank == MPI_WRITER){
                if (iter % PRINTRATE == 0){
                    py::print( "Iter[",iter,"] \t Res = [", buffer_norm , "]\t PhysTime = ", buffer_time);
                }
            }
            // Compute the mean of the norm to stop simulation
            meanNorm = meanNorm/NBEQS;        
        }
    }
    // Write final data
    writeData(iter);
    if (my_rank == MPI_WRITER){
        py::print( "At iter [" , iter ,"] reached norm [" , buffer_norm ,"]\n");
    }

    setResidual(norm_global);
    
    return MeshData::getInstance().get2DData<double>("Cells_py_CC");
}


py::array_t<double> Simulation1D::advance_one(double dt){
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    physTime = 0.;
    
    unsigned int iter = m_timeMethod->getIter();
    py::print("iter = ", iter);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // needed to store the norm and time in a buffer to use py::print
    char buffer_time[100];

    
// iterate till convergence

        // Time discretization
    m_timeMethod->takeStep(dt);

    iter = m_timeMethod->getIter();


    if (iter % SAVERATE == 0 || iter == 1){  writeData(iter);} //write data


    sprintf(buffer_time, "%0.3e",physTime);

    

    return MeshData::getInstance().get2DData<double>("Cells_py_CC");

}
