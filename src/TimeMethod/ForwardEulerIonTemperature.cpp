//
//  ForwardEulerIonTemperature.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ForwardEulerIonTemperature.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("ForwardEulerIonTemperature", ForwardEulerIonTemperature);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include <pybind11/pybind11.h>
#include <mpi.h>

void ForwardEulerIonTemperature::takeStep(double dt)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs = MeshData::getInstance().getData<double>("rhs");
    vector<double>& source = MeshData::getInstance().getData<double>("source");
    //vector<double>& boundaries = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;



    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }

    //compute rhs L2 norm
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        if(iEq == 1 || iEq == 3 || iEq == 4 ){
            m_norm[iEq] = 0;
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/cells[iCell].dx;
                const double rhsU    = rhs[iEq*NBCELLS + iCell];
                const double S_i     = source[iEq*NBCELLS + iCell];
                cells[iCell].uCC[iEq] = cells[iCell].uCC[iEq] - k*rhsU + S_i*dt; // TODO make it general for iEqs
                const double Residual = S_i*cells[0].dx - rhsU;
                m_norm[iEq] += Residual*Residual;
            }
        }
        
    }
    // Space discretization
    m_spaceMethod->discretize();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        if(iEq == 0 || iEq == 2 ){
            m_norm[iEq] = 0;
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/cells[iCell].dx;
                const double rhsU    = rhs[iEq*NBCELLS + iCell];
                const double S_i     = source[iEq*NBCELLS + iCell];
                cells[iCell].uCC[iEq] = cells[iCell].uCC[iEq] - k*rhsU + S_i*dt; // TODO make it general for iEqs
                const double Residual = S_i*cells[0].dx - rhsU;
                m_norm[iEq] += Residual*Residual;
            }
        }
        
    }


    //double k = m_spaceMethod->getDtOvDx();

    //clock_t end = clock();   // end timing
    //py::print( "Time-step on ",  NBCELLS, " cells took ",
    //          double(end-start)/CLOCKS_PER_SEC, " s \n");
    
    physTime += dt;

}
