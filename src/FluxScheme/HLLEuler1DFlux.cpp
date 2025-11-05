//
//  HLLEuler1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 28/01/20.
//  Copyright © 2020 Alejandro Alvarez Laguna. All rights reserved.
//

#include "HLLEuler1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include "../PhysicalModel/EulerEq1D.hpp"
#include <cmath>
using namespace std;

REGISTER_FLUXSCHEME("HLLEuler1D", HLLEuler1DFlux);

vector<double>& HLLEuler1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    // Downcasting as the model is always EulerEq1D
    EulerEq1D* model = dynamic_cast<EulerEq1D*>(m_pm);
    const double gamma  = 5./3.;//model->getGamma();
    // Left variables
    const double rho_L    = uL[0];
    const double m_L      = uL[1];
    const double e_L      = uL[2];
    const double u_L      = m_L/rho_L;
    const double p_L      = (gamma - 1)*(e_L - 0.5*rho_L*pow(u_L,2));
    const double h_L      = (e_L + p_L)/rho_L;
    // Left flux
    const double MassFlux_L   = m_L;
    const double MomFlux_L    = m_L*m_L/rho_L*(3 - gamma)/2 + (gamma - 1)*e_L;
    const double EnergyFlux_L = e_L*m_L/rho_L*gamma - (gamma - 1)/2*pow(m_L,3)/pow(rho_L,2);
    
    // Right variables
    const double rho_R    = uR[0];
    const double m_R      = uR[1];
    const double e_R      = uR[2];
    const double u_R      = m_R/rho_R;
    const double p_R      = (gamma - 1)*(e_R - 0.5*rho_R*pow(u_R,2));
    const double h_R      = (e_R + p_R)/rho_R;
    // Right flux
    const double MassFlux_R   = m_R;
    const double MomFlux_R    = m_R*m_R/rho_R*(3 - gamma)/2 + (gamma - 1)*e_R;
    const double EnergyFlux_R = e_R*m_R/rho_R*gamma - (gamma - 1)/2*pow(m_R,3)/pow(rho_R,2);
    
    // Roe averages
    const double vBar           = (sqrt(rho_L)*u_L + sqrt(rho_R)*u_R)/(sqrt(rho_L) + sqrt(rho_R));
    const double hBar           = (sqrt(rho_L)*h_L + sqrt(rho_R)*h_R)/(sqrt(rho_L) + sqrt(rho_R));
    // Eigenvalues
    const double cBar           = sqrt((gamma - 1)*(hBar - 0.5*pow(vBar,2))); // speed of sound // 1.;
    
//    const double lambda1        = vBar - cBar;
//    const double lambda2        = vBar;
//    const double lambda3        = vBar + cBar;
//    
//    // Eigenvectors
//    const double r1_1 = 1.;
//    const double r1_2 = vBar - cBar;
//    const double r1_3 = hBar -  vBar*cBar;
//    
//    const double r2_1 = 1.;
//    const double r2_2 = vBar;
//    const double r2_3 = vBar*vBar/2.;
//    
//    const double r3_1 = 1.;
//    const double r3_2 = vBar + cBar;
//    const double r3_3 = hBar +  vBar*cBar;
//    
//    const double l1_1 = vBar/(4*cBar)*(2 + (gamma -1)*vBar/cBar);
//    const double l1_2 = -1/(2*cBar)*(1 + (gamma - 1)*vBar/cBar);
//    const double l1_3 = (gamma - 1.)/2.*1/pow(cBar,2);
//    
//    const double l2_1 = 1 - (gamma - 1.)/2.*pow(vBar,2.)/pow(cBar,2.);
//    const double l2_2 = (gamma - 1.)*vBar/pow(cBar,2.);
//    const double l2_3 = -(gamma - 1.)/pow(cBar,2.);
//    
//    const double l3_1 = -vBar/(4*cBar)*(2 - (gamma -1)*vBar/cBar);
//    const double l3_2 = 1/(2*cBar)*(1 - (gamma - 1)*vBar/cBar);
//    const double l3_3 = (gamma - 1.)/2.*1/pow(cBar,2);
//    
//    // Delta values
//    const double deltaMass   = rho_R - rho_L;
//    const double deltaMom    = m_R - m_L;
//    const double deltaEnergy = e_R - e_L;
//    
//    // Numerical dissipation
//    const double alpha1 = l1_1*deltaMass + l1_2*deltaMom + l1_3*deltaEnergy;
//    const double alpha2 = l2_1*deltaMass + l2_2*deltaMom + l2_3*deltaEnergy;
//    const double alpha3 = l3_1*deltaMass + l3_2*deltaMom + l3_3*deltaEnergy;
//    
//    
//    const double numDiss_Mass   = 0.5*(abs(lambda1)*alpha1*r1_1 + abs(lambda2)*alpha2*r2_1 + abs(lambda3)*alpha3*r3_1);
//    const double numDiss_Mom    = 0.5*(abs(lambda1)*alpha1*r1_2 + abs(lambda2)*alpha2*r2_2 + abs(lambda3)*alpha3*r3_2);
//    const double numDiss_Energy = 0.5*(abs(lambda1)*alpha1*r1_3 + abs(lambda2)*alpha2*r2_3 + abs(lambda3)*alpha3*r3_3);
//    
//    // Roe Flux
//    m_flux[0] = 0.5*(MassFlux_R + MassFlux_L) - numDiss_Mass;
//    m_flux[1] = 0.5*(MomFlux_R + MomFlux_L) - numDiss_Mom;
//    m_flux[2] = 0.5*(EnergyFlux_R + EnergyFlux_L) - numDiss_Energy;
    
    const double lambdaBar1        = vBar - cBar;
    const double lambdaBar2        = vBar + cBar;
    const double minLambdaBar      = min(lambdaBar1, lambdaBar2);
    const double maxLambdaBar      = max(lambdaBar1, lambdaBar2);
    
    const double lambdaL1        = u_L - cBar;
    const double lambdaL2        = u_L + cBar;
    const double minLambdaL      = min(lambdaL1, lambdaL2);
    const double maxLambdaL      = max(lambdaL1, lambdaL2);
    
    const double lambdaR1        = u_R - cBar;
    const double lambdaR2        = u_R + cBar;
    const double minLambdaR      = min(lambdaR1, lambdaR2);
    const double maxLambdaR      = max(lambdaR1, lambdaR2);
    
    const double S_L = min(minLambdaBar, min(minLambdaL, minLambdaR));
    const double S_R = max(maxLambdaBar, max(maxLambdaL, maxLambdaR));
    //const double S_L = min(minLambdaL, minLambdaR);
    //const double S_R = max(maxLambdaL, maxLambdaR);
    
    if(S_L > 0)
    {
        m_flux[0] = MassFlux_L;
        m_flux[1] = MomFlux_L;
        m_flux[2] = EnergyFlux_L;
    }
    else if(S_L <= 0 && S_R >= 0)
    {
        m_flux[0] = (S_R*MassFlux_L   - S_L*MassFlux_R   + S_R*S_L*(rho_R - rho_L))/(S_R - S_L);
        m_flux[1] = (S_R*MomFlux_L    - S_L*MomFlux_R    + S_R*S_L*(m_R - m_L))/(S_R - S_L);
        m_flux[2] = (S_R*EnergyFlux_L - S_L*EnergyFlux_R + S_R*S_L*(e_R - e_L))/(S_R - S_L);
    }
    else if(S_R < 0)
    {
        m_flux[0] = MassFlux_R;
        m_flux[1] = MomFlux_R;
        m_flux[2] = EnergyFlux_R;
    }
    // We hardcode the Lax-Friedrich
//    const double maxLambda = max(abs(minLambdaBar), abs(maxLambdaBar));
//    m_flux[0] = 0.5*(MassFlux_R + MassFlux_L) - 0.5*maxLambda*(rho_R - rho_L);
//    m_flux[1] = 0.5*(MomFlux_R + MomFlux_L) - 0.5*maxLambda*(m_R - m_L);
//    m_flux[2] = 0.5*(EnergyFlux_R + EnergyFlux_L) - 0.5*maxLambda*(e_R - e_L);

    return m_flux;
}

