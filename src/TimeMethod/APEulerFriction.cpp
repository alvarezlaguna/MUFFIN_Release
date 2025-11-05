//
//  APEulerFriction.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "APEulerFriction.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("APEulerFriction", APEulerFriction);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include <pybind11/pybind11.h>
#include <mpi.h>

void APEulerFriction::takeStep(double dt)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs = MeshData::getInstance().getData<double>("rhs");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;
    

    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }
    
    //compute rhs L2 norm
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS ; ++iFluid){
        const int NBEQS_per_Fluid = 2;
        const int RHO_i = iFluid*NBEQS_per_Fluid;
        const int RHOU_i = iFluid*NBEQS_per_Fluid + 1;
        for(unsigned int iEq = iFluid*NBEQS_per_Fluid; iEq < (iFluid + 1)*NBEQS_per_Fluid; ++iEq){
            m_norm[iEq] = 0;
            // Euler Friction Parameters
            const double epsilon = EPSILON;
            const double c_sound = SOUNDSPEED[iFluid];
            const double sigma   = 1.;
            for (int iCell = 0; iCell < NBCELLS; ++iCell) {
                double k             = dt/cells[iCell].dx;
                double dx            = cells[iCell].dx;
                const double rhsU    = rhs[iEq*NBCELLS + iCell];        //APEulerFriction1DFlux
                const double S_i     = source[iEq*NBCELLS + iCell];     //EulerFriction1DSourceTerm
                
                // Compute the coefficients for the AP scheme (We assume Neumann conditions)
                double rho_iP1, rho_i, rho_iM1;
                double u_iP1, u_i, u_iM1;
                rho_i   = cells[iCell].uCC[RHO_i];
                rho_iP1 = (iCell == NBCELLS - 1) ? rho_i : cells[iCell + 1].uCC[RHO_i];
                rho_iM1 = (iCell == 0) ?           rho_i : cells[iCell - 1].uCC[RHO_i];
                u_i     = cells[iCell].uCC[RHOU_i]/rho_i;
                u_iP1   = (iCell == NBCELLS - 1) ? u_i : cells[iCell + 1].uCC[RHOU_i]/rho_iP1;
                u_iM1   = (iCell == 0) ?           u_i : cells[iCell - 1].uCC[RHOU_i]/rho_iM1;
                
                const double u_bar_iP12 = (sqrt(rho_i)*u_i + sqrt(rho_iP1)*u_iP1)/(sqrt(rho_i) + sqrt(rho_iP1));
                const double u_bar_iM12 = (sqrt(rho_i)*u_i + sqrt(rho_iM1)*u_iM1)/(sqrt(rho_i) + sqrt(rho_iM1));
                const double b_P_iP12 = u_bar_iP12 + c_sound;
                const double b_P_iM12 = u_bar_iM12 + c_sound;
                const double b_M_iP12 = u_bar_iP12 - c_sound;
                const double b_M_iM12 = u_bar_iM12 - c_sound;
                
                const double alpha   = 2*epsilon*c_sound/(2*epsilon*c_sound + sigma*dx);
                const double beta_j  = (b_P_iM12 - b_M_iP12)/(2*c_sound)*alpha;
                const double gamma_j = 1./(1. + beta_j*sigma*dt/pow(epsilon,2.));
                if(iEq == RHO_i){ // Density
                    cells[iCell].uCC[RHO_i] = cells[iCell].uCC[RHO_i] - alpha*k*rhsU + S_i*dt; // TODO make it general for iEqs
                }
                if(iEq == RHOU_i){ // Momentum
                    cells[iCell].uCC[RHOU_i] = cells[iCell].uCC[RHOU_i] - gamma_j*alpha*k*rhsU + gamma_j*beta_j*S_i*dt; // TODO make it general for iEqs
                }
                const double Residual = S_i*cells[0].dx - rhsU;
                m_norm[iEq] += Residual*Residual;
//                cout<<"uCC["<<iCell<<"].["<<iEq<<"] = "<<cells[iCell].uCC[iEq]<<"\n";
//                cout<<"rhsU["<<iCell<<"].["<<iEq<<"] = "<<rhsU<<"\n";
//                cout<<"S_i["<<iCell<<"].["<<iEq<<"] = "<<S_i<<"\n";
            }
        }
    }
    //clock_t end = clock();   // end timing
    //py::print( "Time-step on ",  NBCELLS, " cells took ",
    //          double(end-start)/CLOCKS_PER_SEC, " s \n");
    
    physTime += dt;

}
