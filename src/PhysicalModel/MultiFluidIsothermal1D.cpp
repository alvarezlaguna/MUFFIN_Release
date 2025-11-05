//
//  MultiFluidIsothermal1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "MultiFluidIsothermal1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("MultiFluidIsothermal1D", MultiFluidIsothermal1D);
#include<algorithm> // for copy() and assign()
#include<iterator> // for back_inserter
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& MultiFluidIsothermal1D::getEigenvalues(const CellDataRef u)
{
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEulerIsoth = 2;
        
        PhysicalModel& pm = *m_physicalModels[iFluid];
        m_uModels[iFluid] = u.withOffset(iFluid*sizeEulerIsoth);
        m_eigenModels[iFluid] = pm.getEigenvalues(m_uModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < sizeEulerIsoth; iEq++){m_eigenvals[iEq + iFluid*sizeEulerIsoth] = m_eigenModels[iFluid][iEq];}
    }
    return m_eigenvals;
}

double MultiFluidIsothermal1D::getMaxEigenvalue(const CellDataRef u)
{
    //cout<<"Entering Max Eigenval\n";
    m_eigenvals = MultiFluidIsothermal1D::getEigenvalues(u);
    
    //cout<<"Max eigenVal = "<<abs(*max_element(m_eigenvals.begin(), m_eigenvals.end()))<<"\n";
    
    return abs(*max_element(m_eigenvals.begin(), m_eigenvals.end()));
}

vector<double>& MultiFluidIsothermal1D::computePhysicalFlux(const CellDataRef u)
{
    
    for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
        unsigned int sizeEulerIsoth = 2;
        
        PhysicalModel& pm = *m_physicalModels[iFluid];
        m_uModels[iFluid] = u.withOffset(iFluid*sizeEulerIsoth);
        m_fluxModels[iFluid] = pm.computePhysicalFlux(m_uModels[iFluid]);
        
        for (unsigned int iEq = 0 ; iEq < m_fluxModels[iFluid].size(); iEq++){m_physFlux[iEq + iFluid*sizeEulerIsoth] = m_fluxModels[iFluid][iEq];}
    }
    return m_physFlux;
}
