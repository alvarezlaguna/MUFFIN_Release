//
//  FourierHeatFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FourierHeatFlux_hpp
#define FourierHeatFlux_hpp

#include <stdio.h>
#include "DiffusiveFlux.hpp"
#include "../CollisionalData/SingleIonMixture.hpp"
#include "../LinearSolver/LinearSolver.hpp"
#include "../LinearSolver/ThomasAlgorithm.hpp"

using namespace Parameters;

class FourierHeatFlux : public DiffusiveFlux
{
public:
    FourierHeatFlux(string name) : DiffusiveFlux(name), m_A(NBCELLS,vector<double>(NBCELLS,0.)), m_B(NBCELLS, 0.), m_x(NBCELLS, 0.) {}
    ~FourierHeatFlux(){}
    
    void setup(){
        DiffusiveFlux::setup();
    
    // copy the boundaries ref
vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");
        m_uInlet  = boundaries[0];
        m_uOutlet = boundaries[1];
    
        m_mixture = static_cast<SingleIonMixture*> (this->getMixture());
        m_mi   = m_mixture->getIonMass();
        m_me   = m_mixture->getElecMass();
        m_ngas = m_mixture->getGasDensity();
        
        m_linearSolver.reset(new ThomasAlgorithm("ThomasAlgorithm", NBCELLS));
        m_linearSolver->setup();

    }
    
    void computeDiffusiveStep(const double dt);
    
    double getConductivity(const CellDataRef u)
    {
        const double e           = 1.602176634e-19;
        const double rho_i       = u[0];
        const double ionMom_z    = u[1];
        const double TotalEnergy = u[2];
        
        const double n      = rho_i/m_mi;
        const double u_xi   = ionMom_z/(rho_i);
        const double p_e    = 2./3.*(TotalEnergy - 0.5*m_mi*n*u_xi*u_xi);
        const double T_e    = p_e/(e*n);
        
        const double Omega22_ee = m_mixture->computeOmega22_ee(T_e);
        const double Omega11_en = m_mixture->computeOmega11_en(T_e);
        const double Omega13_en = m_mixture->computeOmega13_en(T_e);
        
        double nu_eg = m_ngas*(32/15*Omega13_en - 16/3*Omega11_en);
        double nu_ee = n*2./3.*16/5*Omega22_ee;
        
        double conductivity = 5./2.*e/m_me*p_e*1./(nu_eg + nu_ee);
        
        return conductivity;
    }
    
protected:
    SingleIonMixture* m_mixture;
    unique_ptr<LinearSolver> m_linearSolver;
    double m_mi;
    double m_me;
    double m_ngas;
    Matrix m_A;         // Matrix of the linear solver
    vector<double> m_x; // Solution of the linear solver
    vector<double> m_B; // rhs of the linear solver
    CellDataRef m_uInlet;
    CellDataRef m_uOutlet;
    
};

#endif /* FourierHeatFlux_hpp */
