//
//  RoeEulerIsothermal1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "RoeEulerIsothermal1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
using namespace std;

REGISTER_FLUXSCHEME("RoeEulerIsothermal1D", RoeEulerIsothermal1DFlux);

vector<double>& RoeEulerIsothermal1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    // Downcasting as the physical model is always Euler Isothermal
    EulerIsothermal1D *model = dynamic_cast<EulerIsothermal1D*>(m_pm);
    const double soundSpeed  = model->getSoundSpeed();
    
    // Left variables
    const double rho_L    = uL[0];
    const double m_L      = uL[1];
    const double u_L      = m_L/rho_L;
    // Left flux
    const double MassFlux_L   = m_L;
    const double MomFlux_L    = m_L*m_L/rho_L + rho_L*soundSpeed*soundSpeed;
    
    // Right variables
    const double rho_R    = uR[0];
    const double m_R      = uR[1];
    const double u_R      = m_R/rho_R;
    // Right flux
    const double MassFlux_R   = m_R;
    const double MomFlux_R    = m_R*m_R/rho_R + rho_R*soundSpeed*soundSpeed;
    
    // Roe averages
    const double vBar           = (sqrt(rho_L)*u_L + sqrt(rho_R)*u_R)/(sqrt(rho_L) + sqrt(rho_R));
    // Eigenvalues
    const double cBar           = soundSpeed; // speed of sound
    const double lambda1        = vBar - cBar;
    const double lambda2        = vBar + cBar;
    
    // Eigenvectors
    const double r1_1 = 1.;
    const double r1_2 = vBar - cBar;
    
    const double r2_1 = 1.;
    const double r2_2 = vBar + cBar;
    
    const double l1_1 = 1/(2.*soundSpeed)*(vBar + soundSpeed);
    const double l1_2 = -1/(2.*soundSpeed);
    
    const double l2_1 = 1/(2.*soundSpeed)*(-vBar + soundSpeed);
    const double l2_2 = 1/(2.*soundSpeed);
    
    // Delta values
    const double deltaMass   = rho_R - rho_L;
    const double deltaMom    = m_R - m_L;
    
    // Numerical dissipation
    const double alpha1 = l1_1*deltaMass + l1_2*deltaMom;
    const double alpha2 = l2_1*deltaMass + l2_2*deltaMom;
    
    double absLambda1 = abs(lambda1);
    double absLambda2 = abs(lambda2);    
    
    const double numDiss_Mass   = 0.5*(absLambda1*alpha1*r1_1 + absLambda2*alpha2*r2_1);
    const double numDiss_Mom    = 0.5*(absLambda1*alpha1*r1_2 + absLambda2*alpha2*r2_2);
    
    // Roe Flux
    m_flux[0] = 0.5*(MassFlux_R + MassFlux_L) - numDiss_Mass;
    m_flux[1] = 0.5*(MomFlux_R + MomFlux_L) - numDiss_Mom;
    
    return m_flux;
}
