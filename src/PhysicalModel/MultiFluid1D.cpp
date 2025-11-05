//
//  MultiFluid1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 25/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "MultiFluid1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("MultiFluid1D", MultiFluid1D);
#include<algorithm> // for copy() and assign()
#include<iterator> // for back_inserter
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& MultiFluid1D::getEigenvalues(const CellDataRef u)
{
    //cout<<"Entering getEigenvalues\n";
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEuler = 3;
        
        PhysicalModel& pm = *m_physicalModels[iFluid];
        m_uModels[iFluid] = u.withOffset(iFluid*sizeEuler);
        
        m_eigenModels[iFluid] = pm.getEigenvalues(m_uModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < sizeEuler; iEq++){m_eigenvals[iEq + iFluid*sizeEuler] = m_eigenModels[iFluid][iEq];}
    }
    
    return m_eigenvals;
}

double MultiFluid1D::getMaxEigenvalue(const CellDataRef u)
{
    m_eigenvals = MultiFluid1D::getEigenvalues(u);
        
    return abs(*max_element(m_eigenvals.begin(), m_eigenvals.end()));
}

vector<double>& MultiFluid1D::computePhysicalFlux(const CellDataRef u)
{
    
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEuler = 3;
        
        PhysicalModel& pm = *m_physicalModels[iFluid];
        m_uModels[iFluid] = u.withOffset(iFluid*sizeEuler);

        m_fluxModels[iFluid] = pm.computePhysicalFlux(m_uModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < m_fluxModels[iFluid].size(); iEq++){
            m_physFlux[iEq + iFluid*sizeEuler] = m_fluxModels[iFluid][iEq];
        }
    }
    return m_physFlux;
}
