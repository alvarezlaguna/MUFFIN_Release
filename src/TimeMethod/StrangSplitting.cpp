//
//  StrangSplitting.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "StrangSplitting.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("StrangSplitting", StrangSplitting);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

void StrangSplitting::takeStep(double dt)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    //// old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;

    DELTAT    = dt;
    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }
    
    // Initialize the norm with zeroes
    std::fill(m_norm.begin(), m_norm.end(), 0.);
    
    // We take half time step for the source term
    double dt_1 = dt/2.0;
    //First step of the RK
    // The fluxes and the source term are computed in the space method discretization already
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            u_0[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + S_i*dt_1;
        }
    }


    // Second step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = S_i*dt_1;
            (*m_cells)[iCell].uCC[iEq] = 0.75*u_1[iEq*NBCELLS + iCell] + 0.25*(*m_cells)[iCell].uCC[iEq] + 0.25*L_i;
        }
    }

    // Third step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = S_i*dt_1;
            (*m_cells)[iCell].uCC[iEq] = 1./3.*u_1[iEq*NBCELLS + iCell] + 2./3.*(*m_cells)[iCell].uCC[iEq] + 2./3.*L_i;
        }
    }

    // Full step of the fluxes
    // First step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeFluxes();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - k*rhsU;
        }
    }

    // Second step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeFluxes();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = - k*rhsU;
            (*m_cells)[iCell].uCC[iEq] = 0.75*u_1[iEq*NBCELLS + iCell] + 0.25*(*m_cells)[iCell].uCC[iEq] + 0.25*L_i;
        }
    }

    // Third step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeFluxes();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = - k*rhsU;
            (*m_cells)[iCell].uCC[iEq] = 1./3.*u_1[iEq*NBCELLS + iCell] + 2./3.*(*m_cells)[iCell].uCC[iEq] + 2./3.*L_i;
            // Compute the norm
            //const double dU = (*m_cells)[iCell].uCC[iEq] - u_1[iEq*NBCELLS + iCell];

            //m_norm[iEq] += dU*dU/(dt*dt);
        }
    }

    // We take half time step for the source term
    dt_1 = dt / 2.0;
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    //First step of the RK
    // The fluxes and the source term are computed in the space method discretization already
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + S_i*dt_1;
        }
    }


    // Second step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = S_i*dt_1;
            (*m_cells)[iCell].uCC[iEq] = 0.75*u_1[iEq*NBCELLS + iCell] + 0.25*(*m_cells)[iCell].uCC[iEq] + 0.25*L_i;
        }
    }

    // Third step of the RK
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = S_i*dt_1;
            (*m_cells)[iCell].uCC[iEq] = 1./3.*u_1[iEq*NBCELLS + iCell] + 2./3.*(*m_cells)[iCell].uCC[iEq] + 2./3.*L_i;
        }
    }
    

    // Loop to compute norm
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double dU = (*m_cells)[iCell].uCC[iEq] - u_0[iEq*NBCELLS + iCell];
            m_norm[iEq] += dU*dU/(dt*dt);
        }
    }
    
//    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
//        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
//            cout<<"(*m_cells)["<<iCell<<"].uCC["<<iEq<<"] = "<< (*m_cells)[iCell].uCC[iEq]<<"\n";
//        }
//    }
//    
//    for(int i = 0; i<rhs.size(); ++i)   {cout<<"(*m_rhs)["<<i<<"]    = "<<(*m_rhs)[i]<<"\n";}
//    for(int i = 0; i<source.size(); ++i){cout<<"(*m_source)["<<i<<"] = "<<(*m_source)[i]<<"\n";}
//
    physTime += dt;
    
}
