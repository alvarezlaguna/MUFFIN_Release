//
//  FractionalTime1OSplitting.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FractionalTime1OSplitting.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("FractionalTime1OSplitting", FractionalTime1OSplitting);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

void FractionalTime1OSplitting::takeStep(double dt)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    //// old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;
    
    // First step
    // Space discretization


    DELTAT    = dt;
    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }
    
    // Initialize the norm with zeroes
    std::fill(m_norm.begin(), m_norm.end(), 0.);
    
    // We advance the fluxes
    //First step 
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){ 
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            double k             = dt/(*m_cells)[iCell].dx; 
            const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
            u_1[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
            (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - k*rhsU;
        }
    }


    // Second step we advance the momentum source terms
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeSource();
    for (unsigned int iFluid = 0; iFluid < NBFLUIDS; ++iFluid){
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            int iEq = Parameters::MOMENTUMINDICES[iMom];
            // py::print("iEq = ", iEq,"iFluid = ", iFluid); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                const double S_i     = (*m_source).at(iEq, iCell);
                (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + S_i*dt;
            }
        }
    }

    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeFluxes();    
    m_spaceMethod->computeSource();
    // Third step we advance the rest of variables (non-momentum equations)
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        // Check if iEq is a momentum index
        bool isMomentumEq = false;
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            if (iEq == static_cast<unsigned int>(Parameters::MOMENTUMINDICES[iMom])){
                isMomentumEq = true;
                break;
            }
        }
        
        // Only advance fluxes for non-momentum equations
        if (!isMomentumEq) {
            // py::print("Advancing non-momentum equation iEq = ", iEq); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/(*m_cells)[iCell].dx; 
                const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
                const double S_i     = (*m_source).at(iEq, iCell);
                (*m_cells)[iCell].uCC[iEq] = u_1[iEq*NBCELLS + iCell] - k*rhsU + S_i*dt;
            }
        }
    }

    // Loop to compute norm
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double dU = (*m_cells)[iCell].uCC[iEq] - u_1[iEq*NBCELLS + iCell];
            m_norm[iEq] += dU*dU/(dt*dt);
        }
    }
    
    physTime += dt;
    
}
