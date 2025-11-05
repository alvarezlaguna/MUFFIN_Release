//
//  RoeMultiFluid1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "RoeMultiFluid1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
using namespace std;

REGISTER_FLUXSCHEME("RoeMultiFluid1D", RoeMultiFluid1DFlux);

vector<double>& RoeMultiFluid1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEuler = 3;
        
        FluxScheme& fs = *m_fluxSchemes[iFluid];
            m_uRModels[iFluid] = uR.withOffset(iFluid*sizeEuler);
        m_uLModels[iFluid] = uL.withOffset(iFluid*sizeEuler);
        
        m_fluxModels[iFluid] = fs(m_uLModels[iFluid], m_uRModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < m_fluxModels[iFluid].size(); iEq++){
            m_flux[iEq + iFluid*sizeEuler] = m_fluxModels[iFluid][iEq];
        }
    }
    
    return m_flux;
}
