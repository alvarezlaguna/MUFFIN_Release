//
//  SixMomentAQMOM1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "SixMomentAQMOM1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
#include <iostream>

using namespace std;

REGISTER_FLUXSCHEME("SixMomentAQMOM1DFlux", SixMomentAQMOM1DFlux);

vector<double>& SixMomentAQMOM1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    
    double maxEigenVal = max(abs(m_pm->getMaxEigenvalue(uL)),abs(m_pm->getMaxEigenvalue(uR)));
    m_flux_l      = m_pm->computePhysicalFlux(uL);
    // m_eigenvals_l = m_pm->getEigenvalues(uL);
    m_flux_r      = m_pm->computePhysicalFlux(uR);
    // m_eigenvals_r = m_pm->getEigenvalues(uR);

    for (unsigned int iEq = 0; iEq < 5; iEq++){ //1DQMOM System with Lax-Friedrichs flux
        
        m_flux[iEq] = 0.5*(m_flux_l[iEq] + m_flux_r[iEq]) - maxEigenVal*0.5*(uR[iEq] - uL[iEq]);
    }


    // if (m_flux_l[5] > 0 && m_flux_r[5] > 0){
    //     m_flux[5] = m_flux_l[5];
    // }
    // else if (m_flux_r[5] < 0 && m_flux_l[5] < 0){
    //     m_flux[5] = m_flux_r[5];
    // }
    // else{
    //     const double rhoL = uL[0];
    //     const double rhoR = uR[0];  
    //     const double vL = uL[1]/rhoL;
    //     const double vR = uR[1]/rhoR;
    //     const double max_v = max(abs(vL),abs(vR));
    //     m_flux[5] = 0.5*(m_flux_l[5] + m_flux_r[5]) - max_v*0.5*(uR[5] - uL[5]);
    // }
    const double rhoL = uL[0];
    const double rhoR = uR[0];  
    const double vL = uL[1]/rhoL;
    const double vR = uR[1]/rhoR;
    const double max_v = max(abs(vL),abs(vR));
    m_flux[5] = 0.5*(m_flux_l[5] + m_flux_r[5]) - max_v*0.5*(uR[5] - uL[5]);
    // cout<<"Flux = "<<m_flux[5]<<"\n";
    // cout<<"FluxL = "<<m_flux_l[5]<<"\n";
    // cout<<"FluxR = "<<m_flux_r[5]<<"\n";



    return m_flux;

}


