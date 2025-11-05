//
//  SSPRK3.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "SSPRK3.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("SSPRK3", SSPRK3);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

void SSPRK3::takeStep(double dt)
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
    
    //We take the diffusive step
    m_diffFlux->computeDiffusiveStep(dt/2.);
    
    /// STEP WITH JUST THE FLUXES
    STRANG_SPLITTING_STEP = 0;
    
    //First step of the RK
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - 0.5*k*rhsU + 0.5*S_i*dt;
        }
    }


    // Second step of the RK
    m_spaceMethod->discretize();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = - k*rhsU + S_i*dt;
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
        }
    }

    // Third step of the RK
    m_spaceMethod->discretize();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = - k*rhsU + S_i*dt;
            (*m_cells)[iCell].uCC[iEq] = 2./3.*u_1[iEq*NBCELLS + iCell] + 1./3.*(*m_cells)[iCell].uCC[iEq] + 1./6.*L_i;
        }
    }

    // Fourth step of the RK
    m_spaceMethod->discretize();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            const double S_i     = (*m_source).at(iEq, iCell);
            const double L_i     = - k*rhsU  + S_i*dt ;
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
        }
    }

    // /// STEP WITH JUST THE FLUXES
    // STRANG_SPLITTING_STEP = 0;
    
    // //First step of the RK
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - 0.5*k*rhsU + 0.5*S_i*(dt/2.);
    //     }
    // }


    // // Second step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*(dt/2.);
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
    //     }
    // }

    // // Third step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*(dt/2.);
    //         (*m_cells)[iCell].uCC[iEq] = 2./3.*u_1[iEq*NBCELLS + iCell] + 1./3.*(*m_cells)[iCell].uCC[iEq] + 1./6.*L_i;
    //     }
    // }

    // // Fourth step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU  + S_i*(dt/2.) ;
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
    //     }
    // }

    // /// STEP WITH JUST THE SOURCE TERMS
    // STRANG_SPLITTING_STEP = 1;

    // // Space discretization
    // m_spaceMethod->discretize();

    // //First step of the RK
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - 0.5*k*rhsU + 0.5*S_i*dt;
    //     }
    // }


    // // Second step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*dt;
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
    //     }
    // }

    // // Third step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*dt;
    //         (*m_cells)[iCell].uCC[iEq] = 2./3.*u_1[iEq*NBCELLS + iCell] + 1./3.*(*m_cells)[iCell].uCC[iEq] + 1./6.*L_i;
    //     }
    // }

    // /// STEP WITH JUST THE FLUXES AGAIN
    // STRANG_SPLITTING_STEP = 2;

    // // Space discretization
    // m_spaceMethod->discretize();

    // //First step of the RK
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - 0.5*k*rhsU + 0.5*S_i*(dt/2.);
    //     }
    // }


    // // Second step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*(dt/2.);
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
    //     }
    // }

    // // Third step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU + S_i*(dt/2.);
    //         (*m_cells)[iCell].uCC[iEq] = 2./3.*u_1[iEq*NBCELLS + iCell] + 1./3.*(*m_cells)[iCell].uCC[iEq] + 1./6.*L_i;
    //     }
    // }

    // // Fourth step of the RK
    // m_spaceMethod->discretize();
    // for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    //     for (int iCell = 0; iCell < NBCELLS; ++iCell) {
    //         const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
    //         const double S_i     = (*m_source).at(iEq, iCell);
    //         const double L_i     = - k*rhsU  + S_i*(dt/2.) ;
    //         (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*L_i;
    //     }
    // }
    
    m_diffFlux->computeDiffusiveStep(dt/2.);
    // Loop to compute norm
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double dU = (*m_cells)[iCell].uCC[iEq] - u_1[iEq*NBCELLS + iCell];
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
