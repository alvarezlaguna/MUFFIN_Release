//
//  TwoFluidIsothermal1DSourceTermMPI.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermal1DSourceTermMPI.hpp"

// Self-register with the factory
#include "SourceTermRegistrar.hpp"
REGISTER_SOURCETERM("TwoFluidIsothermalMPI1D", TwoFluidIsothermal1DSourceTermMPI);
#include <iostream>
#include <math.h>
#include <cmath>

//#include <petscksp.h>
#include "petsc.h"
#include "petscksp.h"
//#include "petscsys.h"

#include <pybind11/pybind11.h>

//#include <mpi4py/mpi4py.h>

#include <mpi.h>

#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../LinearSolver/ThomasAlgorithm.hpp"

using namespace std;
namespace py = pybind11;

void TwoFluidIsothermal1DSourceTermMPI::setup()
{
    //m_linearSolver.reset(new ThomasAlgorithm("ThomasAlgorithm", NBCELLS));
    //m_linearSolver->setup();
    
    // Set-up the data to store Phi
    MeshData& md = MeshData::getInstance();
    md.createData<double>("Phi", NBCELLS);
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    phi = m_x; // Initializing phi to the initial value that is set in the options
    
    
}

void TwoFluidIsothermal1DSourceTermMPI::solveLinearSystem(Mat& m_A,Vec& m_x,Vec& m_B)
{
    PetscErrorCode ierr;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Solve linear system
     */
    ierr = KSPSolve(m_ksp,m_b_vec,m_x_vec);
    
    // Assemble the vector x
    VecAssemblyBegin(m_x_vec);
    VecAssemblyEnd(m_x_vec);
    
}

double TwoFluidIsothermal1DSourceTermMPI::computeIonizationConstant()
{
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");
    vector<double>& x      = MeshData::getInstance().getData<double>("x");

    int size;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double T_e     = 1.;
    const double ionizConst = IONIZCONST;
    const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    double Dx    = cells[0].dx;
    const int n1 = NBCELLS - 1;
    double integral = 0.;
    
    // First cell
    // Order of boundaries
    // boundaries = [rho_eInlet, rho_eOutlet, rhoU_eInlet, rhoU_eOutlet, rho_iInlet, rho_iOutlet, rhoU_iInlet, rhoU_iOutlet,]
    double rhoe_Ghost = boundaries[0][0];
    double rho_iP1    = cells[0].uCC[0];
    double rhoe_Wall  = (rhoe_Ghost + rho_iP1)/2.;
    integral += (rhoe_Wall + rho_iP1)/2.*Dx/2.; // We sum the first half cell
    

    for (unsigned int iCell = 0; iCell < n1; ++iCell)
    {
        Dx    = x[iCell + 1] - x[iCell];
        double rho_i = cells[iCell].uCC[0];
        rho_iP1      = cells[iCell + 1].uCC[0];
        integral += (rho_iP1 + rho_i)/2.*Dx; // We sum the first half cell
    }
    // Last cell
    Dx          = cells[n1].dx;
    rhoe_Ghost  = boundaries[1][0];
    double rho_i       = cells[n1].uCC[0];
    rhoe_Wall  = (rhoe_Ghost + rho_i)/2.;
    integral += (rhoe_Wall + rho_i)/2.*Dx/2.; // We sum the first half cell
    integral = integral*ionization_Constant;
    
    // Sum over the processors
    double totalIntegral = 0;
    MPI_Allreduce(&integral, &totalIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    double Da = 1;
    /// Flux of ions
    if(rank == size - 1)
    {
        double Flux_i = abs(boundaries[1][3]);
        Da = 2*Flux_i/totalIntegral;
    }
    
    MPI_Bcast(&Da, 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
    
    if (Da > 0) {// At the beginning the flux can become negative and then the ionization also is negative. We set it to zero in this case.
        return Da;
    }
    else {
        return 0.;
    }
    
}

void TwoFluidIsothermal1DSourceTermMPI::computeSource()
{
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
    vector<double>& x      = MeshData::getInstance().getData<double>("x");
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    const double Dx = cells[0].dx; //Assuming the cells have the same delta_x
    const double DebyeLength = DEBYELENGTH;
    const double MassRatio   = MASSRATIO;
    // Initializing the frequency
    double frequency = 0.;
    int    rank, size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int nGlobal = 0;
    for(int iP = 0; iP < size; ++iP){
        nGlobal += NBCELLS_MPI[iP];
    }
    
    // Compute the matrix for the Poisson solver

    //PetscPrintf(MPI_COMM_WORLD,"Before Loop\n");
    int offsetRank = 0;
    for (int iP = 0; iP < rank; ++iP){
        offsetRank += NBCELLS_MPI[iP];
    }
    int offset = 0;
    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
    {
        offset = offsetRank + iCell;
        const double rho_e   = cells[iCell].uCC[0];
        const double rho_i   = cells[iCell].uCC[2];
        
        m_B[iCell] = DebyeLength*DebyeLength*(rho_e - rho_i)*Dx*Dx; //
        
        // Insert these values to the PETSc vector
        const double value = m_B[iCell];

        VecSetValues(m_b_vec,1,&offset,&value,INSERT_VALUES);

    }
    //cout<<"rank = "<<rank<<"\t NBCELLS = "<<NBCELLS<<"\t Dx = "<<Dx<<"\n";
    // View the vector
    //VecView(m_b_vec,PETSC_VIEWER_STDOUT_WORLD);
    // Add the boundary values for the Poisson Solver
    const double Phi_in  = 2*PHIIN - phi[0];
    const double Phi_out = 2*PHIOUT - phi[NBCELLS - 1];
    
    // The first rank adds the intlet conditions
    if(rank == 0){
        int firstCell  = 0;
        
        double firstBoundary = -Phi_in;

        VecSetValues(m_b_vec,1,&firstCell,&firstBoundary,ADD_VALUES);

    }
    // The last rank adds the outlet conditions
    if(rank == size - 1){
        int lastCell   = nGlobal - 1;
        
        double lastBoundary  = -Phi_out;
        
        VecSetValues(m_b_vec,1,&lastCell,&lastBoundary,ADD_VALUES);
        
    }
    //PetscPrintf(MPI_COMM_WORLD,"Before Assembly\n");
    // Assemble the vector to use it
    VecAssemblyBegin(m_b_vec);
    VecAssemblyEnd(m_b_vec);

    // Compute the variable Phi (that is in m_x)
    solveLinearSystem(m_A_matrix, m_x_vec, m_b_vec);
    //PetscPrintf(MPI_COMM_WORLD,"After Solve\n");

    
    // We assume that the order of the vars is
    // {rho_e, rho_eu_e, rho_i, rho_iu_i}
    // {0,     1,        2,     3       }
    // Compute ionization
    const double ionizConst = IONIZCONST;
    const double T_e        = 1.;
    const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    const double nn         = computeIonizationConstant();
    const double nu_iz      = nn*ionization_Constant;
    
    offsetRank = 0;
    for (int iP = 0; iP < rank; ++iP){
        offsetRank += NBCELLS_MPI[iP];
    }
    offset = 0;
    
    double phiZeroMOne = 0; // Inlet ghost cell
    double phiZero = 0;     // First cell in processor
    double phiLastPOne = 0; // Outlet ghost cell
    double phiLast = 0;     // Last cell in processor
    
    int left_neigh = (rank>0) ? rank-1 : MPI_PROC_NULL;
    int right_neigh = (rank<(size-1)) ? rank+1 : MPI_PROC_NULL;
    MPI_Status status;
    
    // Get and send the first cell to the other processor
    VecGetValues(m_x_vec,1,&offsetRank,&phiZero);
    MPI_Sendrecv(&phiZero, 1, MPI_DOUBLE, left_neigh, 6, &phiLastPOne,  1, MPI_DOUBLE, right_neigh, 6, MPI_COMM_WORLD, &status);
    // Get and send the last cell to the other processor
    int offsetLast = offsetRank + NBCELLS - 1;
    VecGetValues(m_x_vec,1,&offsetLast,&phiLast);
    MPI_Sendrecv(&phiLast, 1, MPI_DOUBLE, right_neigh, 6, &phiZeroMOne,  1, MPI_DOUBLE, left_neigh, 6, MPI_COMM_WORLD, &status);
    // Put ghost cells in first and last processor
    if(rank == 0) phiZeroMOne = Phi_in;
    if(rank == size - 1) phiLastPOne = Phi_out;
    
    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
    {
        offset = offsetRank + iCell;

        // Compute the variables
        const double rho_e   = cells[iCell].uCC[0];
        const double rho_i   = cells[iCell].uCC[2];
        //gradPhi = mdVec(1).gradient(1,:);
        const double u_e     = cells[iCell].uCC[1]/rho_e;
        const double u_i     = cells[iCell].uCC[3]/rho_i;
        const double massRatio = MASSRATIO;
        const double collConstIons = COLLIONS;
        const double collConstElec = COLLELECTRONS;
        
        // Compute Ionization Constant with an integration
        const double ionization = nu_iz*rho_e;
        
        // compute grad Phi
//        double gradPhi = 0.;
//        double denominator, firstTerm, secondTerm, thirdTerm;
//        if(iCell == 0) // First cell
//        {
//            x_im1       = x[0] - cells[0].dx;
//            x_i         = x[0];
//            x_ip1       = x[1];
//            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
//
//            firstTerm  = -Phi_in*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
//            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
//            thirdTerm  = m_x[iCell + 1]*(x_i - x_im1)*(x_i - x_im1)/denominator;
//
//            gradPhi = firstTerm + secondTerm + thirdTerm;
//
//        }
//        else if (iCell == NBCELLS - 1) // Last cell
//        {
//            x_im1       = x[iCell - 1];
//            x_i         = x[iCell];
//            x_ip1       = x[n1] + cells[n1].dx;
//            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
//
//            firstTerm  = -m_x[iCell - 1]*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
//            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
//            thirdTerm  = Phi_out*(x_i - x_im1)*(x_i - x_im1)/denominator;
//
//            gradPhi = firstTerm + secondTerm + thirdTerm;
//
//        }
//        else { // Inner cells
//            x_im1       = x[iCell - 1];
//            x_i         = x[iCell];
//            x_ip1       = x[iCell + 1];
//            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
//
//            firstTerm = -m_x[iCell - 1]*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
//            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
//            thirdTerm = m_x[iCell + 1]*(x_i - x_im1)*(x_i - x_im1)/denominator;
//
//            gradPhi = firstTerm + secondTerm + thirdTerm;
//        }
        // compute grad Phi
        double gradPhi = 0.;
        VecGetValues(m_x_vec,1,&offset,&phi[iCell]);
        double phi_nMone;
        double phi_nPone;

        if(iCell == 0) // First cell
        {
            //get the value at iCell = 1
            int offsetPone = offset + 1;
            VecGetValues(m_x_vec,1,&offsetPone,&phi_nPone);
            gradPhi = (phi_nPone - phiZeroMOne)/(2*Dx);

        }
        else if (iCell == NBCELLS - 1){ // Last cell {
            //get the value at iCell = i - 1
            int offsetMone = offset - 1;
            VecGetValues(m_x_vec,1,&offsetMone,&phi_nMone);

            gradPhi = (phiLastPOne - phi_nMone)/(2*Dx);
        }
        else {
            //get the value at iCell = i - 1
            int offsetMone = offset - 1;
            VecGetValues(m_x_vec,1,&offsetMone,&phi_nMone);
            //get the value at iCell = i + 1
            int offsetPone = offset + 1;
            VecGetValues(m_x_vec,1,&offsetPone,&phi_nPone);
            gradPhi = (phi_nPone - phi_nMone)/(2*Dx);
        }
        source[0*NBCELLS + iCell] = ionization;
        source[1*NBCELLS + iCell] = (gradPhi*rho_e*massRatio - collConstElec*rho_e*u_e);
        source[2*NBCELLS + iCell] = ionization;
        source[3*NBCELLS + iCell] = (-gradPhi*rho_i - collConstIons*rho_i*u_i);
        
        // compute frequency
        const double freqIzElec     = ionization/rho_e;
        const double freqCollElec   = collConstElec;
        const double freqCollIons   = collConstIons;
        const double freqPlasma     = sqrt(rho_e)*DebyeLength*sqrt(MassRatio);
        const double maxFreq_iCell  = max(freqIzElec,max(freqCollElec,max(freqCollIons,freqPlasma)));
        frequency = max(frequency,maxFreq_iCell);
    }
    
    double totalFrequency = 0;
    MPI_Allreduce(&frequency, &totalFrequency, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //PetscPrintf(MPI_COMM_WORLD,"MPI_Allreduce\n");
    
    // Set the Frequency in the base class
    SourceTerm::setFrequency(frequency);
    //for(int i = 0; i < source.size(); ++i){cout<<"source = "<<source[i]<<"\n";}
    
}
