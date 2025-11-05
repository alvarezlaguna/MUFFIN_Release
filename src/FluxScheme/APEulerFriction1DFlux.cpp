//
//  APEulerFriction1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "APEulerFriction1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"

REGISTER_FLUXSCHEME("APEulerFriction1D", APEulerFriction1DFlux);

vector<double>& APEulerFriction1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
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
    
    // HLL averages
    const double vBar           = (sqrt(rho_L)*u_L + sqrt(rho_R)*u_R)/(sqrt(rho_L) + sqrt(rho_R));
    // Eigenvalues
    const double cBar              = soundSpeed; // speed of sound
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
    
    //const double S_L = min(minLambdaBar, min(minLambdaL, minLambdaR));
    //const double S_R = max(maxLambdaBar, max(maxLambdaL, maxLambdaR));
    const double S_L = minLambdaBar;
    const double S_R = maxLambdaBar;
    
    if(S_L > 0)
    {
        m_flux[0] = MassFlux_L;
        m_flux[1] = MomFlux_L;
    }
    else if(S_L <= 0 && S_R >= 0)
    {
        m_flux[0] = (S_R*MassFlux_R - S_L*MassFlux_L - (MomFlux_R - MomFlux_L))/(S_R - S_L);
        m_flux[1] = (S_R*MomFlux_L - S_L*MomFlux_R + S_R*S_L*(m_R - m_L))/(S_R - S_L);
    }
    else if(S_R < 0)
    {
        m_flux[0] = MassFlux_R;
        m_flux[1] = MomFlux_R;
    }

    m_flux[0] *= 1/EPSILON;
    m_flux[1] *= 1/EPSILON;
    
    return m_flux;
}
