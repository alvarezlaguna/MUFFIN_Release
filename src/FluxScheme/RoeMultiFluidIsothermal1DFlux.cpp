//
//  RoeMultiFluidIsothermal1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 02/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "RoeMultiFluidIsothermal1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
using namespace std;

REGISTER_FLUXSCHEME("RoeMultiFluidIsothermal1D", RoeMultiFluidIsothermal1DFlux);

vector<double>& RoeMultiFluidIsothermal1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEulerIsoth = 2;
        
        FluxScheme& fs = *m_fluxSchemes[iFluid];
        m_uRModels[iFluid] = uR.withOffset(iFluid*sizeEulerIsoth);
        m_uLModels[iFluid] = uL.withOffset(iFluid*sizeEulerIsoth);
        
        m_fluxModels[iFluid] = fs(m_uLModels[iFluid], m_uRModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < m_fluxModels[iFluid].size(); iEq++){
            m_flux[iEq + iFluid*sizeEulerIsoth] = m_fluxModels[iFluid][iEq];
        }
    }
    
    return m_flux;
}
