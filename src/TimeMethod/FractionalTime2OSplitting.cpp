//
//  FractionalTime2OSplitting.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FractionalTime2OSplitting.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("FractionalTime2OSplitting", FractionalTime2OSplitting);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <cmath>
#include <mpi.h>

void FractionalTime2OSplitting::setup()
{
    // Call parent class setup
    TimeMethod::setup();
    
    // Setup EFieldPhi source if it exists
    if (EFIELDFRACTIONAL && MeshData::getInstance().hasData("EFieldPhi")) {
        m_EFieldSource = &MeshData::getInstance().get2DData<double>("EFieldPhi");
    } else {
        m_EFieldSource = nullptr;
    }
}

void FractionalTime2OSplitting::takeStep(double dt)
{

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    //// old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;
    
    // First step
    // Space discretization it's done 


    DELTAT    = dt;
    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }
    
    // Initialize the norm with zeroes
    std::fill(m_norm.begin(), m_norm.end(), 0.);
    
    // We advance the fluxes
    // Copy all cell data to u_0 before advancing
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            u_0[iEq*NBCELLS + iCell] = (*m_cells)[iCell].uCC[iEq];
        }
    }
    // copy of EField source at time n
    std::memcpy(E_0.data(), m_EFieldSource->data(0, 0), NBCELLS * sizeof(double));
    // copy of Phi source at time n
    std::memcpy(Phi_0.data(), m_EFieldSource->data(1, 0), NBCELLS * sizeof(double));

    //First step just advance the densities
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        // Check if iEq is a density index (momentum - 1)
        bool isDensityEq = false;
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            if (iEq == static_cast<unsigned int>(Parameters::MOMENTUMINDICES[iMom] - 1)){
                isDensityEq = true;
                break;
            }
        }
        // Only advance fluxes for density equations
        if (isDensityEq) {
            // py::print("Advancing density equation iEq = ", iEq); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/(*m_cells)[iCell].dx; 
                const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
                (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] - 0.5*k*rhsU;
            }
        }
    }

    // Second step we compute the EField at time n+1/2 and advance the momentum equations
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeEField();
    // Update the EField source to time n+1/2
    // Get pointers once for all EField and Phi operations
    double* EField_ptr = m_EFieldSource->mutable_data(0, 0);
    double* Phi_ptr = m_EFieldSource->mutable_data(1, 0);
    for (int iCell = 0; iCell < NBCELLS; ++iCell) {
        double E_half = 0.5 * (E_0[iCell] + EField_ptr[iCell]);
        EField_ptr[iCell] = E_half;
        double Phi_half = 0.5 * (Phi_0[iCell] + Phi_ptr[iCell]);
        Phi_ptr[iCell] = Phi_half;
    }

    //Third step we advance the momentum equations with this EField at time n+1/2
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        // Check if iEq is a density index (momentum - 1)
        bool isDensityEq = false;
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            if (iEq == static_cast<unsigned int>(Parameters::MOMENTUMINDICES[iMom] - 1)){
                isDensityEq = true;
                break;
            }
        }
        
        // Only advance source terms for non-density equations
        if (!isDensityEq) {
            // py::print("Advancing non-density equation iEq = ", iEq); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/(*m_cells)[iCell].dx; 
                const double S_i = (*m_source).at(iEq, iCell);
                const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];

                (*m_cells)[iCell].uCC[iEq] = (*m_cells)[iCell].uCC[iEq] + 0.5*S_i*dt - 0.5*k*rhsU;
            }
        }
    }

    // Fourth step compute the fluxes and just advance again the densities
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeFluxes();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        // Check if iEq is a density index (momentum - 1)
        bool isDensityEq = false;
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            if (iEq == static_cast<unsigned int>(Parameters::MOMENTUMINDICES[iMom] - 1)){
                isDensityEq = true;
                break;
            }
        }
        // Only advance fluxes for density equations
        if (isDensityEq) {
            // py::print("Advancing density equation AGAIN iEq = ", iEq); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/(*m_cells)[iCell].dx; 
                const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];
                (*m_cells)[iCell].uCC[iEq] = u_0[iEq*NBCELLS + iCell] - k*rhsU;
            }
        }
    }

    // Fifth step we compute the EField at time n+1 and advance the momentum equations
    // copy of EField source at time n+1/2
    // std::memcpy(E_0.data(), m_EFieldSource->data(0, 0), NBCELLS * sizeof(double));
    m_spaceMethod->setBoundaries();
    m_spaceMethod->computeEField();
    // Update the EField source to time n+1 (reuse the same pointers)
    for (int iCell = 0; iCell < NBCELLS; ++iCell) {
        double E_half = 0.5 * (E_0[iCell] + EField_ptr[iCell]);
        EField_ptr[iCell] = E_half;
        double Phi_half = 0.5 * (Phi_0[iCell] + Phi_ptr[iCell]);
        Phi_ptr[iCell] = Phi_half;
    }

    //Third step we advance the momentum equations with this EField at time n+1
    m_spaceMethod->computeSource();
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        // Check if iEq is a density index (momentum - 1)
        bool isDensityEq = false;
        for (unsigned int iMom = 0; iMom < Parameters::MOMENTUMINDICES.size(); ++iMom){
            if (iEq == static_cast<unsigned int>(Parameters::MOMENTUMINDICES[iMom] - 1)){
                isDensityEq = true;
                break;
            }
        }
        
        // Only advance source terms for non-density equations
        if (!isDensityEq) {
            // py::print("Advancing non-density equation AGAIN iEq = ", iEq); //--- IGNORE ---
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/(*m_cells)[iCell].dx; 
                const double S_i     = (*m_source).at(iEq, iCell);
                const double rhsU    = (*m_rhs)[iEq*NBCELLS + iCell];

                (*m_cells)[iCell].uCC[iEq] = u_0[iEq*NBCELLS + iCell] + S_i*dt - k*rhsU;
            }
        }
    }
    
    // Loop to compute norm
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double dU = (*m_cells)[iCell].uCC[iEq] - u_0[iEq*NBCELLS + iCell];
            m_norm[iEq] += dU*dU/(dt*dt);
        }
    }
    
    physTime += dt;
    
    
}
