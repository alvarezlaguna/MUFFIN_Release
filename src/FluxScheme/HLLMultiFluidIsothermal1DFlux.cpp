//
//  HLLMultiFluidIsothermal1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/01/20.
//  Copyright © 2020 Alejandro Alvarez Laguna. All rights reserved.
//

#include "HLLMultiFluidIsothermal1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
using namespace std;

REGISTER_FLUXSCHEME("HLLMultiFluidIsothermal1D", HLLMultiFluidIsothermal1DFlux);

vector<double>& HLLMultiFluidIsothermal1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEulerIsoth = 2;
        
        FluxScheme& fs = *m_fluxSchemes[iFluid];
        m_uRModels[iFluid] = uR.withOffset(iFluid*sizeEulerIsoth);
        m_uLModels[iFluid] = uL.withOffset(iFluid*sizeEulerIsoth);
        
        m_fluxModels[iFluid] = fs(m_uLModels[iFluid], m_uRModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < sizeEulerIsoth; iEq++){
            m_flux[iEq + iFluid*sizeEulerIsoth] = m_fluxModels[iFluid][iEq];
        }
    }
    
    return m_flux;
}
