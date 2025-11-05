//
//  TwoFluidIsothermal1DSourceTerm.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermal1DSourceTerm.hpp"

// Self-register with the factory
#include "SourceTermRegistrar.hpp"
REGISTER_SOURCETERM("TwoFluidIsothermal1D", TwoFluidIsothermal1DSourceTerm);
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

void TwoFluidIsothermal1DSourceTerm::setup()
{
    m_linearSolver.reset(new ThomasAlgorithm("ThomasAlgorithm", NBCELLS));
    m_linearSolver->setup();
    
    // Set-up the data to store Phi
    MeshData& md = MeshData::getInstance();
    md.createData<double>("Phi", NBCELLS);
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    phi = m_x; // Initializing phi to the initial value that is set in the options
    
    
}

double TwoFluidIsothermal1DSourceTerm::computeIonizationConstant()
{
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");

    vector<double>& x      = MeshData::getInstance().getData<double>("x");
    
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
    
    // Flux of ions
    double Flux_i = abs(boundaries[1][3]);
    double Da = 2*Flux_i/integral;
    
    // Flux of the electrons
    //double soundSpeed_e = SOUNDSPEED[0]; //We assume the electrons to be the first fluid
    //double rhoe_eG      = cells[0].uCC[0];
    //const double pi     = atan(1)*4;
    //double FluxW_e      = rhoe_eG*soundSpeed_e/sqrt(2*pi);
    //double Da = 2*FluxW_e/integral;
    
    if (Da > 0) {// At the beginning the flux can become negative and then the ionization also is negative. We set it to zero in this case.
        return Da;
    }
    else {
        return 0.;
    }
    
}

void TwoFluidIsothermal1DSourceTerm::computeSource()
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
    
    //m_linearSolver->solveLinearSystem(m_A, m_x, m_B);
    // Compute the matrix for the Poisson solver
    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
    {
        const double rho_e   = cells[iCell].uCC[0];
        const double rho_i   = cells[iCell].uCC[2];
        
        m_B[iCell] = DebyeLength*DebyeLength*(rho_e - rho_i);
    }
    // Add the boundary values for the Poisson Solver
    const double Phi_in = 2*PHIIN - phi[0];
    const double Phi_out = 2*PHIOUT - phi[NBCELLS -1];
    unsigned int n1 = NBCELLS -1;
    // First cell
    double x_im1    = x[0] - cells[0].dx;
    double x_i      = x[0];
    double x_ip1    = x[1];
    double denominator = 1./2.*(x_ip1 - x_im1)*(x_ip1 - x_i)*(x_i - x_im1);

    m_B[0] = m_B[0] - Phi_in*(x_ip1 - x_i)/denominator;
    
    // Last Cell
    x_im1    = x[n1 - 1];
    x_i      = x[n1];
    x_ip1    = x[n1] + cells[n1].dx;
    denominator = 1./2.*(x_ip1 - x_im1)*(x_ip1 - x_i)*(x_i - x_im1);
    m_B[NBCELLS - 1] = m_B[NBCELLS - 1] - Phi_out*(x_i - x_im1)/denominator;
    
//    for (int iCell = 0; iCell< NBCELLS; iCell++){
//        cout<<m_B[iCell]<<"\n";
//    }
    
    // Compute the variable Phi (that is in m_x)
    m_linearSolver->solveLinearSystem(m_A, m_x, m_B);
    phi = m_x; // copy the value of m_x to phi

    
    // We assume that the order of the vars is
    // {rho_e, rho_eu_e, rho_i, rho_iu_i}
    // {0,     1,        2,     3       }
    // Compute ionization
    const double ionizConst = IONIZCONST;
    const double T_e        = 1.;
    const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    const double nn         = computeIonizationConstant();
    const double nu_iz      = nn*ionization_Constant;
    
    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
    {
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
        double gradPhi = 0.;
        double denominator, firstTerm, secondTerm, thirdTerm;
        if(iCell == 0) // First cell
        {
            x_im1       = x[0] - cells[0].dx;
            x_i         = x[0];
            x_ip1       = x[1];
            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
            
            firstTerm  = -Phi_in*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
            thirdTerm  = m_x[iCell + 1]*(x_i - x_im1)*(x_i - x_im1)/denominator;
            
            gradPhi = firstTerm + secondTerm + thirdTerm;

        }
        else if (iCell == NBCELLS - 1) // Last cell
        {
            x_im1       = x[iCell - 1];
            x_i         = x[iCell];
            x_ip1       = x[n1] + cells[n1].dx;
            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
            
            firstTerm  = -m_x[iCell - 1]*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
            thirdTerm  = Phi_out*(x_i - x_im1)*(x_i - x_im1)/denominator;
            
            gradPhi = firstTerm + secondTerm + thirdTerm;

        }
        else { // Inner cells
            x_im1       = x[iCell - 1];
            x_i         = x[iCell];
            x_ip1       = x[iCell + 1];
            denominator = (x_i - x_im1)*(x_ip1 - x_i)*(x_ip1 - x_im1);
            
            firstTerm = -m_x[iCell - 1]*(x_ip1 - x_i)*(x_ip1 - x_i)/denominator;
            secondTerm = m_x[iCell]*(x_ip1 + x_im1 - 2*x_i)*(x_ip1 - x_im1)/denominator;
            thirdTerm = m_x[iCell + 1]*(x_i - x_im1)*(x_i - x_im1)/denominator;
            
            gradPhi = firstTerm + secondTerm + thirdTerm;
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
    // Set the Frequency in the base class
    SourceTerm::setFrequency(frequency);
    //for(int i = 0; i < source.size(); ++i){cout<<"source = "<<source[i]<<"\n";}
    
}
