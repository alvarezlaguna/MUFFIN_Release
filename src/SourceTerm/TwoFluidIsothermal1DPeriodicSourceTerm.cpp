//
//  TwoFluidIsothermal1DPeriodicSourceTerm.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 06/02/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermal1DPeriodicSourceTerm.hpp"

// Self-register with the factory
#include "SourceTermRegistrar.hpp"
REGISTER_SOURCETERM("TwoFluidIsothermal1DPeriodic", TwoFluidIsothermal1DPeriodicSourceTerm);
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
#include "../LinearSolver/ThomasAlgorithmPeriodic.hpp"

using namespace std;
namespace py = pybind11;

void TwoFluidIsothermal1DPeriodicSourceTerm::setup()
{
    m_linearSolver.reset(new ThomasAlgorithmPeriodic("ThomasAlgorithmPeriodic", NBCELLS));
    m_linearSolver->setup();
    
    // Set-up the data to store Phi
    MeshData& md = MeshData::getInstance();
    md.createData<double>("Phi", NBCELLS);
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    phi = m_x; // Initializing phi to the initial value that is set in the options
    
    
}

double TwoFluidIsothermal1DPeriodicSourceTerm::computeIonizationConstant()
{
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");

    const double T_e     = 1.;
    const double ionizConst = IONIZCONST;
    const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    const double Dx = cells[0].dx; //Assuming the cells have the same delta_x
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
        double rho_i = cells[iCell].uCC[0];
        rho_iP1      = cells[iCell + 1].uCC[0];
        integral += (rho_iP1 + rho_i)/2.*Dx; // We sum the first half cell
    }
    // Last cell
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

void TwoFluidIsothermal1DPeriodicSourceTerm::computeSource()
{
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
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
        
        m_B[iCell] = DebyeLength*DebyeLength*(rho_e - rho_i)*Dx*Dx; // The
    }
    // We set the last value to zero
    m_B[NBCELLS-1] = 0;
    
    // Compute the variable Phi (that is in m_x)
    m_linearSolver->solveLinearSystem(m_A, m_x, m_B);
    phi = m_x; // copy the value of m_x to phi
    
    
    // We assume that the order of the vars is
    // {rho_e, rho_eu_e, rho_i, rho_iu_i}
    // {0,     1,        2,     3       }
    // Compute ionization
    //const double ionizConst = IONIZCONST;
    //const double T_e        = 1.;
    //const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    //const double nn         = computeIonizationConstant();
    //const double nu_iz      = nn*ionization_Constant;
    
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
        //const double ionization = nu_iz*rho_e;
        
        // compute grad Phi
        double gradPhi = 0.;
        if(iCell == 0) // First cell
        {
            gradPhi = (m_x[1] - m_x[NBCELLS - 1])/(2*Dx);
            //cout<<"gradPhi[0] = "<<gradPhi<<"\n";
        }
        else if (iCell == NBCELLS - 1) // Last cell
        {
            gradPhi = (m_x[0] - m_x[NBCELLS - 2])/(2*Dx);
            //cout<<"gradPhi[n1] = "<<gradPhi<<"\n";
        }
        else { // Inned cells
            gradPhi = (m_x[iCell + 1] - m_x[iCell - 1])/(2*Dx);
        }
        source[0*NBCELLS + iCell] = 0;
        source[1*NBCELLS + iCell] = (gradPhi*rho_e*massRatio - collConstElec*rho_e*u_e);
        source[2*NBCELLS + iCell] = 0;
        source[3*NBCELLS + iCell] = (-gradPhi*rho_i - collConstIons*rho_i*u_i);
        
        // compute frequency
        //const double freqIzElec     = ionization/rho_e;
        const double freqCollElec   = collConstElec;
        const double freqCollIons   = collConstIons;
        const double freqPlasma     = sqrt(rho_e)*DebyeLength*sqrt(MassRatio);
        const double maxFreq_iCell  = max(freqCollElec,max(freqCollIons,freqPlasma));
        frequency = max(frequency,maxFreq_iCell);
    }
    // Set the Frequency in the base class
    SourceTerm::setFrequency(frequency);
    //for(int i = 0; i < source.size(); ++i){cout<<"source = "<<source[i]<<"\n";}
    
}
