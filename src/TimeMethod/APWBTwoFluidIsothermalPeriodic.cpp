//
//  APWBTwoFluidIsothermalPeriodic.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 30/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "APWBTwoFluidIsothermalPeriodic.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("APWBTwoFluidIsothermalPeriodic", APWBTwoFluidIsothermalPeriodic);
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../PhysicalModel/PhysicalModelFactory.hpp"
#include "../SourceTerm/SourceTermFactory.hpp"
#include "../BoundaryCondition/BoundaryConditionFactory.hpp"
#include "../SpaceReconstructor/ReconstructorFactory.hpp"
#include "../FluxScheme/FluxSchemeFactory.hpp"
#include "../LinearSolver/ThomasAlgorithmPeriodic.hpp"
#include <iostream>
#include <cmath>

void APWBTwoFluidIsothermalPeriodic::setup(){
    
    py::print("Setting-up:\t ",this->getName(),"\n");
    
    // TODO: Change this
    vector<CellDataRef>& bounds = MeshData::getInstance().getData<CellDataRef>("boundaries");
    m_uInlet  = bounds[0];
    m_uOutlet = bounds[1];
    
    // Create the Physical model
    m_pm = PhysicalModelFactory::CreatePhysicalModel(PHYSICALMODELNAME);
    
    // Initialize the flux scheme
    m_flux = FluxSchemeFactory::CreateFluxScheme(FLUXSCHEMENAME);
    m_flux->setPhysicalModel(m_pm.get());
    
    // Initialize the source term
    m_source = SourceTermFactory::CreateSourceTerm(SOURCETERMNAME);
    m_source->setup();
    m_source->setPhysicalModel(m_pm.get());
    
    // Initialize the reconstructor for second order solution in space
    m_reconstructor = ReconstructorFactory::CreateReconstructor(RECONSTRUCTORNAME);
    m_reconstructor->setup(); // initialize the limiter
    
    // TODO: Change this to be set from input file
    m_InletBC  = BoundaryConditionFactory::CreateBoundaryCondition(INLETTYPE,"Left");
    
    m_OutletBC  = BoundaryConditionFactory::CreateBoundaryCondition(OUTLETTYPE,"Right");
    
    py::print("Space discretization using the flux scheme:\t ",m_flux->getName(),"\n");
    py::print("Space discretization using the source term:\t ",m_source->getName(),"\n");
    py::print("Space discretization using the reconstructor:\t ",m_reconstructor->getName(),"\n");
    
    // Setting up linear solver
    m_linearSolver.reset(new ThomasAlgorithmPeriodic("ThomasAlgorithmPeriodic", NBCELLS));
    m_linearSolver->setup();
    
    // Set-up the data to store Phi
    MeshData& md = MeshData::getInstance();
    md.createData<double>("Phi", NBCELLS);
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    phi = m_x; // Initializing phi to the initial value that is set in the options
    
}

double APWBTwoFluidIsothermalPeriodic::velocityFlux(const CellDataRef uL, const CellDataRef uR) {
    const double n_eL = uL[0];
    const double n_eR = uR[0];
    const double v_eL = uL[1]/n_eL;
    const double v_eR = uR[1]/n_eR;
    
    //const double u_12 = (sqrt(n_eL)*v_eL + sqrt(n_eR)*v_eR)/(sqrt(n_eL) + sqrt(n_eR));
    const double u_12 = (v_eL + v_eR)/2.;
    
    //Standard
    //return 0.5*(v_eL + v_eR) - 0.5*sqrt(MASSRATIO)*(n_eR - n_eL);
    
    //Balanced
    double Mach12 = abs(u_12/sqrt(MASSRATIO));
    double M_co   = (0./sqrt(MASSRATIO));
    double M_o    = min(1., max(Mach12, M_co));
    double fa     = sqrt((1 - M_o)*(1 - M_o)*Mach12*Mach12 + 4*M_o*M_o)/(1 + M_o*M_o);
    
    return 0.5*(v_eL + v_eR) - 0.5*sqrt(MASSRATIO)*fa*(1./n_eR - 1./n_eL)*(n_eR + n_eL)/2.;
    
    // Old implementation
    //const double Mach12 = max(abs(u_12/sqrt(MASSRATIO)), 1.);
    //return 0.5*(v_eL + v_eR) - 0.5/(sqrt(MASSRATIO)*Mach12)*(n_eR - n_eL);
    
}

double APWBTwoFluidIsothermalPeriodic::velocityFlux_Euler(const CellDataRef uL, const CellDataRef uR) {
    const double n_eL = uL[0];
    const double n_eR = uR[0];
    const double v_eL = uL[1]/n_eL;
    const double v_eR = uR[1]/n_eR;
    
    const double u_12 = (sqrt(n_eL)*v_eL + sqrt(n_eR)*v_eR)/(sqrt(n_eL) + sqrt(n_eR));
    
    //return u_12;
    
    double Mach12 = abs(u_12/sqrt(MASSRATIO));
    double M_co   = (0./sqrt(MASSRATIO));
    double M_o    = min(1., max(Mach12, M_co));
    double fa     = sqrt((1 - M_o)*(1 - M_o)*Mach12*Mach12 + 4*M_o*M_o)/(1 + M_o*M_o);
    
    return 0.5*(v_eL + v_eR) - 0.5*sqrt(MASSRATIO)*fa*(n_eR - n_eL);
    
}

double APWBTwoFluidIsothermalPeriodic::pressureFlux(const CellDataRef uL, const CellDataRef uR) {
    
    const double n_eL = uL[0];
    const double n_eR = uR[0];
    const double v_eL = uL[1]/n_eL;
    const double v_eR = uR[1]/n_eR;
    
    const double u_12 = (sqrt(n_eL)*v_eL + sqrt(n_eR)*v_eR)/(sqrt(n_eL) + sqrt(n_eR));
    
    // Standard
    //return 0.5*MASSRATIO*(n_eR + n_eL) - 0.5*sqrt(MASSRATIO)*(v_eR - v_eL);
    // Balanced
    double Mach12 = abs(u_12/sqrt(MASSRATIO));
    double M_co   = (0./sqrt(MASSRATIO));
    double M_o    = min(1., max(Mach12, M_co));
    double fa     = sqrt((1 - M_o)*(1 - M_o)*Mach12*Mach12 + 4*M_o*M_o)/(1 + M_o*M_o);
    
    return 0.5*MASSRATIO*(n_eR + n_eL) - 0.5*sqrt(MASSRATIO)*fa*(v_eR - v_eL)*(n_eR + n_eL)/2.;
    
}

double APWBTwoFluidIsothermalPeriodic::densityNumericalViscosity(const double niP1,const double ni, const double niM1){
    const double pressureViscosity = std::log(niP1) + std::log(niM1) - 2*std::log(ni);
    //const double pressureViscosity = 2*(niP1 - ni)/(niP1 + ni) - 2*(ni - niM1)/(ni + niM1);
    const double Dx = getDx();
    const double LorentzViscosity  = Dx*Dx*DEBYELENGTH*DEBYELENGTH*ni;
    
    const double DtovDx = getDtOvDx();
    
    return DtovDx*DtovDx*MASSRATIO*(pressureViscosity - LorentzViscosity);
}

void APWBTwoFluidIsothermalPeriodic::computeElectronDensity()
{
    setBoundaries();
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    // Set the Dx so it is easier to access later
    const double Dx     = cells[0].dx; //Assuming the cells have the same delta_x
    setDx(Dx);
    const double DtOvDx = getDtOvDx();
    const double Dt     = DtOvDx*Dx;
    
    

    
    // Reconstruction
    m_reconstructor->reconstructField();
    
    
    // Main loop
    
    double n_ii, n_ei;
    double n_eiP1, n_eiM1;
    double u_ip12, u_im12;
    double numVisc, factor;
    
    // First cell
    n_ii = cells[0].uCC[2];
    n_ei = cells[0].uCC[0];
    n_eiP1 = cells[1].uCC[0];
    n_eiM1 = m_uInlet[0];

    u_ip12   = velocityFlux(cells[0].uR, cells[1].uL);
    u_im12   = velocityFlux(m_uInlet, cells[0].uL);
    numVisc     = densityNumericalViscosity(n_eiP1, n_ei, n_eiM1);
    const double theta = 1;
    factor      = (1 + theta*Dt*Dt*n_ii*MASSRATIO*DEBYELENGTH*DEBYELENGTH);
    
    m_inverseDensity[0] = (1/n_ei)*(1 + DtOvDx*(u_ip12 - u_im12) - theta*numVisc)/factor;
    
    int n1 = NBCELLS - 1;
    
    for (unsigned int iCell = 1; iCell < n1; ++iCell)
    {
        n_ii = cells[iCell].uCC[2];
        n_ei = cells[iCell].uCC[0];
        n_eiP1 = cells[iCell + 1].uCC[0];
        n_eiM1 = cells[iCell - 1].uCC[0];
        u_ip12   = velocityFlux(cells[iCell].uR, cells[iCell + 1].uL);
        u_im12   = velocityFlux(cells[iCell - 1].uR, cells[iCell].uL);
        numVisc     = densityNumericalViscosity(n_eiP1, n_ei, n_eiM1);
        factor      = (1 + theta*Dt*Dt*n_ii*MASSRATIO*DEBYELENGTH*DEBYELENGTH);
        
        m_inverseDensity[iCell] = (1/n_ei)*(1 + DtOvDx*(u_ip12 - u_im12) - theta*numVisc)/factor;
        
    }
    
    // Last cell
    n_ii = cells[n1].uCC[2];
    n_ei = cells[n1].uCC[0];
    n_eiP1 = m_uOutlet[0];
    n_eiM1 = cells[n1 - 1].uCC[0];
    u_ip12   = velocityFlux(cells[n1].uR, m_uOutlet);
    u_im12   = velocityFlux(cells[n1 - 1].uR, cells[n1].uL);
    numVisc     = densityNumericalViscosity(n_eiP1, n_ei, n_eiM1);
    factor      = (1 + theta*Dt*Dt*n_ii*MASSRATIO*DEBYELENGTH*DEBYELENGTH);
    
    m_inverseDensity[n1] = (1/n_ei)*(1 + DtOvDx*(u_ip12 - u_im12) - theta*numVisc)/factor;
    
    for (int iCell = 0; iCell<NBCELLS; iCell++)
    {
        // Update of the electron density
        m_density_nP1Minus[iCell] = 1/m_inverseDensity[iCell];
//        cout<<"m_density_nP1Minus["<<iCell<<"] = "<<m_density_nP1Minus[iCell]<<"\n";
    }
    
}

void APWBTwoFluidIsothermalPeriodic::computeElectricPotential(){
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    const double Dx = cells[0].dx; //Assuming the cells have the same delta_x
    const double DebyeLength = DEBYELENGTH;
    
    // Compute the matrix for the Poisson solver
    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
    {
        //const double rho_e   = cells[iCell].uCC[0];
        const double rho_e   = m_density_nP1Minus[iCell];
        const double rho_i   = cells[iCell].uCC[2];
        
        m_B[iCell] = DebyeLength*DebyeLength*(rho_e - rho_i)*Dx*Dx;
    }
    // We set the last value to zero
    //m_B[NBCELLS-1] = 0;
    
    // Compute the variable Phi (that is in m_x)
    m_linearSolver->solveLinearSystem(m_A, m_x, m_B);
    phi = m_x; // copy the value of m_x to phi
    
//    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
//    {
//        cout<<"phi["<<iCell<<"] = "<<m_x[iCell]<<"\n";
//    }

    
    
    // compute grad Phi
    const int n1 = NBCELLS - 1;
    m_gradPhi[0] = (m_x[1] - m_x[NBCELLS - 1])/(2*Dx);
    for(unsigned int iCell = 1; iCell < n1; ++iCell){
        m_gradPhi[iCell] = (m_x[iCell + 1] - m_x[iCell - 1])/(2*Dx);
    }
    
    m_gradPhi[n1] = (m_x[0] - m_x[NBCELLS - 2])/(2*Dx);
}

void APWBTwoFluidIsothermalPeriodic::computeElectronVelocity(){
    
    setBoundaries();
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    
    // TODO: See if I can store the CC values in a reference vector
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    // Set the Dx so it is easier to access later
    const double Dx     = cells[0].dx; //Assuming the cells have the same delta_x
    setDx(Dx);
    const double DtOvDx = getDtOvDx();
    const double Dt     = DtOvDx*Dx;
    

    
    // Reconstruction
    // WATCH OUT!!!!!!!!!!!!!
    // This can be done more efficiently if I reconstruct only the electron density
    m_reconstructor->reconstructField();
    
    // Main loop
    
    double u_ei, n_ei;
    double p_ip12, p_im12;
    
    // First cell
    n_ei = cells[0].uCC[0];
    u_ei = cells[0].uCC[1]/n_ei;
    p_ip12   = pressureFlux(cells[0].uR, cells[1].uL);
    p_im12   = pressureFlux(m_uInlet, cells[0].uL);
    
    m_velocity[0] = u_ei - DtOvDx/n_ei*(p_ip12 - p_im12) + Dt*MASSRATIO*m_gradPhi[0];
    m_momentum_nP1Minus[0] = m_density_nP1Minus[0]*m_velocity[0];
    
    int n1 = NBCELLS - 1;
    
    for (unsigned int iCell = 1; iCell < n1; ++iCell)
    {
        n_ei = cells[iCell].uCC[0];
        u_ei = cells[iCell].uCC[1]/n_ei;
        p_ip12   = pressureFlux(cells[iCell].uR, cells[iCell + 1].uL);
        p_im12   = pressureFlux(cells[iCell - 1].uR, cells[iCell].uL);
        
        m_velocity[iCell] = u_ei - DtOvDx/n_ei*(p_ip12 - p_im12) + Dt*MASSRATIO*m_gradPhi[iCell];
        m_momentum_nP1Minus[iCell] = m_density_nP1Minus[iCell]*m_velocity[iCell];
    }
    
    // Last cell
    n_ei = cells[n1].uCC[0];
    u_ei = cells[n1].uCC[1]/n_ei;
    p_ip12   = pressureFlux(cells[n1].uR, m_uOutlet);
    p_im12   = pressureFlux(cells[n1 - 1].uR, cells[n1].uL);
    
    m_velocity[n1] = u_ei - DtOvDx/n_ei*(p_ip12 - p_im12) + Dt*MASSRATIO*m_gradPhi[n1];
    m_momentum_nP1Minus[n1] = m_density_nP1Minus[n1]*m_velocity[n1];
    
//    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
//    {
//        cout<<"m_velocity["<<iCell<<"] = "<<m_velocity[iCell]<<"\n";
//        cout<<"m_momentum_nP1Minus["<<iCell<<"] = "<<m_momentum_nP1Minus[iCell]<<"\n";
//    }
    
    // Correction step to ensure mass Conservation.
    // First cell
    double u_ip12   = velocityFlux_Euler(cells[0].uR, cells[1].uL);
    double u_im12   = velocityFlux_Euler(m_uInlet, cells[0].uL);
    m_density_nP1Minus[0] = m_density_nP1Minus[0]/(1 + DtOvDx*(u_ip12 - u_im12));

    for (unsigned int iCell = 1; iCell < n1; ++iCell)
    {
        u_ip12   = velocityFlux_Euler(cells[iCell].uR, cells[iCell + 1].uL);
        u_im12   = velocityFlux_Euler(cells[iCell - 1].uR, cells[iCell].uL);
        m_density_nP1Minus[iCell] = m_density_nP1Minus[iCell]/(1 + DtOvDx*(u_ip12 - u_im12));
    }
    // Last cell
    u_ip12   = velocityFlux_Euler(cells[n1].uR, m_uOutlet);
    u_im12   = velocityFlux_Euler(cells[n1 - 1].uR, cells[n1].uL);
    m_density_nP1Minus[n1] = m_density_nP1Minus[n1]/(1 + DtOvDx*(u_ip12 - u_im12));
    
}

void APWBTwoFluidIsothermalPeriodic::computeEulerianStep(){
    
    setBoundaries();
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    
    // TODO: See if I can store the CC values in a reference vector
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    // Set the Dx so it is easier to access later
    const double Dx     = cells[0].dx; //Assuming the cells have the same delta_x
    setDx(Dx);
    const double DtOvDx = getDtOvDx();
    //const double Dt     = DtOvDx*Dx;
    

    // Reconstruction
    // WATCH OUT!!!!!!!!!!!!!
    //m_reconstructor->reconstructField();
    
    // Main loop
    
    double u_ei, n_ei, nu_ei;
    double n_eip12, n_eim12;
    double n_eiP1, n_eiM1;
    double nu_eiP1, nu_eiM1;
    double nu_eip12, nu_eim12;
    double u_ip12, u_im12;

    int n1 = NBCELLS - 1;
    
    // First cell
    n_ei    = m_density_nP1Minus[0];
    n_eiP1  = m_density_nP1Minus[1];
    n_eiM1  = m_density_nP1Minus[n1];
    nu_ei   = m_momentum_nP1Minus[0];
    nu_eiP1 = m_momentum_nP1Minus[1];
    nu_eiM1 = m_momentum_nP1Minus[n1];
    u_ei = nu_ei/n_ei;
    u_ip12   = velocityFlux_Euler(cells[0].uR, cells[1].uL);
    u_im12   = velocityFlux_Euler(m_uInlet, cells[0].uL);
    
    // We store in the previous vectors
    n_eip12 = (u_ip12 > 0)? n_ei: n_eiP1;
    n_eim12 = (u_im12 > 0)? n_eiM1: n_ei;
    nu_eip12 = (u_ip12 > 0)? nu_ei: nu_eiP1;
    nu_eim12 = (u_im12 > 0)? nu_eiM1: nu_ei;
    
    // Step taking only np1_minus
    m_density_nP1[0]  = m_density_nP1Minus[0] - DtOvDx*(n_eip12*u_ip12 - n_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_density_nP1Minus[0];
    m_momentum_nP1[0] = m_momentum_nP1Minus[0] - DtOvDx*(nu_eip12*u_ip12 - nu_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_momentum_nP1Minus[0];
    
    for (unsigned int iCell = 1; iCell < n1; ++iCell)
    {
        n_ei = m_density_nP1Minus[iCell];
        n_eiP1 = m_density_nP1Minus[iCell + 1];
        n_eiM1 = m_density_nP1Minus[iCell - 1];
        nu_eiP1 = m_momentum_nP1Minus[iCell + 1];
        nu_eiM1 = m_momentum_nP1Minus[iCell - 1];
        nu_ei = m_momentum_nP1Minus[iCell];
        u_ei    = nu_ei/n_ei;
        u_ip12   = velocityFlux_Euler(cells[iCell].uR, cells[iCell + 1].uL);
        u_im12   = velocityFlux_Euler(cells[iCell - 1].uR, cells[iCell].uL);
        
        // We store in the previous vectors
        n_eip12 = (u_ip12 > 0)? n_ei: n_eiP1;
        n_eim12 = (u_im12 > 0)? n_eiM1: n_ei;
        nu_eip12 = (u_ip12 > 0)? nu_ei: nu_eiP1;
        nu_eim12 = (u_im12 > 0)? nu_eiM1: nu_ei;

        // Step taking only np1_minus
        m_density_nP1[iCell]  = m_density_nP1Minus[iCell] - DtOvDx*(n_eip12*u_ip12 - n_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_density_nP1Minus[iCell];
        m_momentum_nP1[iCell] = m_momentum_nP1Minus[iCell] - DtOvDx*(nu_eip12*u_ip12 - nu_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_momentum_nP1Minus[iCell];
    }
    
    // Last cell
    n_ei    = m_density_nP1Minus[n1];
    n_eiP1  = m_density_nP1Minus[0];
    n_eiM1  = m_density_nP1Minus[n1 - 1];
    nu_ei   = m_momentum_nP1Minus[n1];
    nu_eiP1 = m_momentum_nP1Minus[0];
    nu_eiM1 = m_momentum_nP1Minus[n1 - 1];
    u_ei     = nu_ei/n_ei;
    u_ip12   = velocityFlux_Euler(cells[n1].uR, m_uOutlet);
    u_im12   = velocityFlux_Euler(cells[n1 - 1].uR, cells[n1].uL);
    
    // We store in the previous vectors
    n_eip12 = (u_ip12 > 0)? n_ei: n_eiP1;
    n_eim12 = (u_im12 > 0)? n_eiM1: n_ei;
    nu_eip12 = (u_ip12 > 0)? nu_ei: nu_eiP1;
    nu_eim12 = (u_im12 > 0)? nu_eiM1: nu_ei;

    m_density_nP1[n1]  = m_density_nP1Minus[n1] - DtOvDx*(n_eip12*u_ip12 - n_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_density_nP1Minus[n1];
    m_momentum_nP1[n1] = m_momentum_nP1Minus[n1] - DtOvDx*(nu_eip12*u_ip12 - nu_eim12*u_im12) + DtOvDx*(u_ip12 - u_im12)*m_momentum_nP1Minus[n1];
    
//    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
//    {
//        cout<<"m_density_nP1["<<iCell<<"] = "<<m_density_nP1[iCell]<<"\n";
//        cout<<"m_momentum_nP1["<<iCell<<"] = "<<m_momentum_nP1[iCell]<<"\n";
//    }
    
}

double APWBTwoFluidIsothermalPeriodic::computeIonizationConstant(){
    
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
    
    
    if (Da > 0) {// At the beginning the flux can become negative and then the ionization also is negative. We set it to zero in this case.
        return Da;
    }
    else {
        return 0.;
    }
    
}

double APWBTwoFluidIsothermalPeriodic::WBmassSourceL(const double u_12, const double c_12){
    
    double firstPart  = 1/(2*c_12)*signfunction(c_12 + u_12);
    double secondPart = 1/(2*c_12)*signfunction(u_12 - c_12);
    
    return 0.5*(firstPart - secondPart);
    
}

double APWBTwoFluidIsothermalPeriodic::WBmassSourceR(const double u_12, const double c_12){
    
    double firstPart  = 1/(2*c_12)*signfunction(c_12 + u_12);
    double secondPart = 1/(2*c_12)*signfunction(u_12 - c_12);
    
    return 0.5*(-firstPart + secondPart);
    
}

double APWBTwoFluidIsothermalPeriodic::WBmomSourceL(const double u_12, const double c_12){
    
    double firstPart  = u_12/(2*c_12)*(signfunction(c_12 + u_12) - signfunction(u_12 - c_12));
    double secondPart = 1/2*(signfunction(c_12 + u_12) + signfunction(u_12 - c_12));
    
    return 0.5*(1 + firstPart + secondPart);
}

double APWBTwoFluidIsothermalPeriodic::WBmomSourceR(const double u_12, const double c_12){
    
    double firstPart  = u_12/(2*c_12)*(signfunction(c_12 + u_12) - signfunction(u_12 - c_12));
    double secondPart = 1/2*(signfunction(c_12 + u_12) + signfunction(u_12 - c_12));
    
    
    return 0.5*(1 - firstPart - secondPart);
}


void APWBTwoFluidIsothermalPeriodic::computeSources(){
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
    vector<double>& phi    = MeshData::getInstance().getData<double>("Phi");
    
    
    const double Dx = cells[0].dx; //Assuming the cells have the same delta_x

    
    // Compute ionization
    //const double ionizConst = IONIZCONST;
    //const double T_e        = 1.;
    //const double ionization_Constant = ionizConst*exp(-EPSIONIZ/T_e);
    //const double nn         = computeIonizationConstant();
    //const double nu_iz      = nn*ionization_Constant;
    
    // First cell
    // Compute the variables
    double rho_e   = cells[0].uCC[0];
    double rho_i   = cells[0].uCC[2];
    double u_e     = cells[0].uCC[1]/rho_e;
    double u_i     = cells[0].uCC[3]/rho_i;
    double massRatio = MASSRATIO;
    double collConstIons = COLLIONS;
    double collConstElec = COLLELECTRONS;
    
    // Compute the interface values
    double rho_i_im1    = m_uInlet[2];
    double rho_i_i      = cells[0].uCC[2];
    double rho_i_ip1    = cells[1].uCC[2];
    double rho_i_ip12   = (rho_i_ip1 + rho_i_i)/2.;
    double rho_i_im12   = (rho_i_im1 + rho_i_i)/2.;
    double u_i_im1   = m_uInlet[3]/m_uInlet[2];
    double u_i_i     = cells[0].uCC[3]/cells[0].uCC[2];
    double u_i_ip1   = cells[1].uCC[3]/cells[1].uCC[2];
    double u_i_im12  = (sqrt(rho_i_im1)*u_i_im1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_im1) + sqrt(rho_i_i));
    double u_i_ip12  = (sqrt(rho_i_ip1)*u_i_ip1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_ip1) + sqrt(rho_i_i));
    // WATCH OUT: hard coded Dirichtlet boundary and delta x
    double gradPhi_im12  = (phi[0] - phi[NBCELLS - 1])/Dx;
    double gradPhi_ip12  = (phi[1] - phi[0])/Dx;
    double sound_i = SOUNDSPEED[1];
    // WB Lorentz contribution
    double MassLorentz = gradPhi_im12*rho_i_im12*WBmassSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmassSourceR(u_i_ip12, sound_i);
    double MomLorentz  = gradPhi_im12*rho_i_im12*WBmomSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmomSourceR(u_i_ip12, sound_i);

    source[0*NBCELLS + 0] = 0.;
    source[1*NBCELLS + 0] = (- collConstElec*rho_e*u_e);
    source[2*NBCELLS + 0] = - MassLorentz;
    source[3*NBCELLS + 0] = (-MomLorentz - collConstIons*rho_i*u_i);
    
    for (unsigned int iCell = 1; iCell < (NBCELLS - 1); ++iCell)
    {
        // Compute the variables
        rho_e   = cells[iCell].uCC[0];
        rho_i   = cells[iCell].uCC[2];
        u_e     = cells[iCell].uCC[1]/rho_e;
        u_i     = cells[iCell].uCC[3]/rho_i;
        massRatio = MASSRATIO;
        collConstIons = COLLIONS;
        collConstElec = COLLELECTRONS;
        
        // Compute the interface values
        rho_i_im1    = cells[iCell - 1].uCC[2];
        rho_i_i      = cells[iCell].uCC[2];
        rho_i_ip1    = cells[iCell + 1].uCC[2];
        rho_i_ip12   = (rho_i_ip1 + rho_i_i)/2.;
        rho_i_im12   = (rho_i_im1 + rho_i_i)/2.;
        u_i_im1   = cells[iCell - 1].uCC[3]/cells[iCell - 1].uCC[2];
        u_i_i     = cells[iCell].uCC[3]/cells[iCell].uCC[2];
        u_i_ip1   = cells[iCell + 1].uCC[3]/cells[iCell + 1].uCC[2];
        u_i_im12  = (sqrt(rho_i_im1)*u_i_im1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_im1) + sqrt(rho_i_i));
        u_i_ip12  = (sqrt(rho_i_ip1)*u_i_ip1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_ip1) + sqrt(rho_i_i));
        
        gradPhi_im12  = (phi[iCell] - phi[iCell - 1])/Dx;
        gradPhi_ip12  = (phi[iCell + 1] - phi[iCell])/Dx;
        sound_i = SOUNDSPEED[1];
        // WB Lorentz contribution
        MassLorentz = gradPhi_im12*rho_i_im12*WBmassSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmassSourceR(u_i_ip12, sound_i);
        MomLorentz  = gradPhi_im12*rho_i_im12*WBmomSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmomSourceR(u_i_ip12, sound_i);
        
        source[0*NBCELLS + iCell] = 0.;
        source[1*NBCELLS + iCell] = (- collConstElec*rho_e*u_e);
        source[2*NBCELLS + iCell] = 0. - MassLorentz;
        source[3*NBCELLS + iCell] = (-MomLorentz - collConstIons*rho_i*u_i);
    }
    // Last cell
    // Compute the variables
    const int n1 = (NBCELLS - 1);
    rho_e   = cells[n1].uCC[0];
    rho_i   = cells[n1].uCC[2];
    u_e     = cells[n1].uCC[1]/rho_e;
    u_i     = cells[n1].uCC[3]/rho_i;
    massRatio = MASSRATIO;
    collConstIons = COLLIONS;
    collConstElec = COLLELECTRONS;
    
    // Compute the interface values
    rho_i_im1    = cells[n1 - 1].uCC[2];
    rho_i_i      = cells[n1].uCC[2];
    rho_i_ip1    = m_uOutlet[2];
    rho_i_ip12   = (rho_i_ip1 + rho_i_i)/2.;
    rho_i_im12   = (rho_i_im1 + rho_i_i)/2.;
    u_i_im1   = cells[n1 - 1].uCC[3]/cells[n1 - 1].uCC[2];
    u_i_i     = cells[n1].uCC[3]/cells[n1].uCC[2];
    u_i_ip1   = m_uOutlet[3]/m_uOutlet[2];
    u_i_im12  = (sqrt(rho_i_im1)*u_i_im1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_im1) + sqrt(rho_i_i));
    u_i_ip12  = (sqrt(rho_i_ip1)*u_i_ip1 + sqrt(rho_i_i)*u_i_i)/(sqrt(rho_i_ip1) + sqrt(rho_i_i));
    
    gradPhi_im12  = (phi[n1] - phi[n1 - 1])/Dx;
    gradPhi_ip12  = (phi[0] - phi[n1])/Dx;
    sound_i = SOUNDSPEED[1];
    // WB Lorentz contribution
    MassLorentz = gradPhi_im12*rho_i_im12*WBmassSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmassSourceR(u_i_ip12, sound_i);
    MomLorentz  = gradPhi_im12*rho_i_im12*WBmomSourceL(u_i_im12, sound_i) + gradPhi_ip12*rho_i_ip12*WBmomSourceR(u_i_ip12, sound_i);
    
    source[0*NBCELLS + n1] = 0.;
    source[1*NBCELLS + n1] = (- collConstElec*rho_e*u_e);
    source[2*NBCELLS + n1] = 0. - MassLorentz;
    source[3*NBCELLS + n1] = (-MomLorentz - collConstIons*rho_i*u_i);
    
}


void APWBTwoFluidIsothermalPeriodic::setBoundaries(){
    
    // old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    m_uInlet  = m_InletBC->setBoundary();
    m_uOutlet = m_OutletBC->setBoundary();
    // TODO: Write it general in a class

}

void APWBTwoFluidIsothermalPeriodic::computeIonFlux(){
    
    setBoundaries();
    
    const int    n1 = NBCELLS - 1;          // number of cells minus 1
    
    FluxScheme& flux1D          = *m_flux;      // reference is used to keep the same syntax as before
    //SourceTerm& sourceterm      = *m_source;
    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs         = MeshData::getInstance().getData<double>("rhs");
    // old boundary call  = MeshData::getInstance().getData<double>("boundaries");
    
    

    
    // Reconstruction
    m_reconstructor->reconstructField();
    
    //double maxEigenvalue = 0;
    
    // Main loop
    m_Fip12         = flux1D(cells[0].uR, cells[1].uL);
    m_Fim12         = flux1D(m_uInlet, cells[0].uL);
    // First cell
    for (unsigned int iEq = 2; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS]        = m_Fip12[iEq] - m_Fim12[iEq];
    }
    
    // Loop over inner cells
    for (int iCell = 1; iCell < n1; ++iCell) {
        m_Fip12         = flux1D(cells[iCell].uR, cells[iCell+1].uL);
        m_Fim12         = flux1D(cells[iCell-1].uR, cells[iCell].uL);
        for (unsigned int iEq = 2; iEq < NBEQS; iEq++){
            // compute the flux
            rhs[iEq*NBCELLS + iCell]    = m_Fip12[iEq] - m_Fim12[iEq];
        }
        
    } // integral flux on internal cells
    // Last cell
    m_Fip12         = flux1D(cells[n1].uR, m_uOutlet);
    m_Fim12         = flux1D(cells[n1-1].uR, cells[n1].uL);
    for (unsigned int iEq = 2; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS + n1]       = m_Fip12[iEq] - m_Fim12[iEq];
    }
    
//    for (unsigned int iCell = 0; iCell < 4*NBCELLS; ++iCell)
//    {
//        cout<<"rhs["<<iCell<<"] = "<<rhs[iCell]<<"\n";
//    }

}
void APWBTwoFluidIsothermalPeriodic::takeStep(double dt){
    APWBTwoFluidIsothermalPeriodic::takeStep();
}

void APWBTwoFluidIsothermalPeriodic::takeStep()
{
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs = MeshData::getInstance().getData<double>("rhs");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    setDt();
    
    m_iter++;
    
    // First step
    // Compute electron density at t = t^*
    computeElectronDensity();
    // Compute Electric Potential
    computeElectricPotential();
    // Compute Electron Velocity
    computeElectronVelocity();
    // Computing electron eulerian step
    computeEulerianStep();
    // Computing sources
    computeSources();
    // Computing the ion convective terms
    computeIonFlux();
    
    
    double k = getDtOvDx();
    double dt = k*cells[0].dx; //Assuming the cells have the same delta_x
    
    // Initialize the norm with zeroes
    std::fill(m_norm.begin(), m_norm.end(), 0.);
    double Residual_Density_e, Residual_Momentum_e;
    
    // Update the electrons
    for (int iCell = 0; iCell < NBCELLS; ++iCell) {
        Residual_Density_e  = (m_density_nP1[iCell] + dt*source[0*NBCELLS + iCell] - cells[iCell].uCC[0])/k;
        Residual_Momentum_e = (m_momentum_nP1[iCell] + dt*source[1*NBCELLS + iCell] - cells[iCell].uCC[1])/k;
        
        // Update
        cells[iCell].uCC[0] = m_density_nP1[iCell] + dt*source[0*NBCELLS + iCell];
        cells[iCell].uCC[1] = m_momentum_nP1[iCell] + dt*source[1*NBCELLS + iCell];
        
        m_norm[0] += Residual_Density_e*Residual_Density_e;
        m_norm[1] += Residual_Momentum_e*Residual_Momentum_e;
    }
    
//    for (unsigned int iEq = 0; iEq < NBEQS; ++iEq)
//    {
//        cout<<"m_norm["<<iEq<<"] = "<<m_norm[iEq]<<"\n";
//    }
    
//    // Check that the norm is not zero
//    if(m_norm[0] != 0 ){ //check that the residual is not -inf
//        m_norm[0] = log10(sqrt(m_norm[0]));
//    }
//    else{ m_norm[0] = 0.; } //when the residual is 0, we set the log to 0 as well
//    if(m_norm[1] != 0 ){ //check that the residual is not -inf
//        m_norm[1] = log10(sqrt(m_norm[1]));
//    }
//    else{ m_norm[1] = 0.; } //when the residual is 0, we set the log to 0 as well
    
    
    //Update the ions
    for(unsigned int iEq = 2; iEq < NBEQS; ++iEq){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            const double rhsU    = rhs[iEq*NBCELLS + iCell];
            const double S_i     = source[iEq*NBCELLS + iCell];
            
            cells[iCell].uCC[iEq] = cells[iCell].uCC[iEq] - k*rhsU + S_i*dt;
            const double Residual = S_i*cells[0].dx - rhsU;
            m_norm[iEq] += Residual*Residual;
        }
        // Check that the norm is not zero
//        if(m_norm[iEq] != 0 ){ //check that the residual is not -inf
//            m_norm[iEq] = log10(sqrt(m_norm[iEq]));
//        }
//        else{ m_norm[iEq] = 0.; } //when the residual is 0, we set the log to 0 as well
    }
    
//    for (unsigned int iCell = 0; iCell < 4*NBCELLS; ++iCell)
//    {
//        cout<<"rhs["<<iCell<<"] = "<<rhs[iCell]<<"\n";
//    }
//
//    for (unsigned int iEq = 0; iEq < NBEQS; ++iEq)
//    {
//        for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
//        {
//            cout<<"source["<<iEq<<"]["<<iCell <<"] = "<<source[iEq*NBCELLS + iCell]<<"\n";
//        }
//        for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell)
//        {
//            cout<<"cells["<<iEq<<"]["<<iCell <<"] = "<<cells[iCell].uCC[iEq]<<"\n";
//        }
//    }
    
    

    physTime += dt;
    py::print("dt = ",dt);
    
}
