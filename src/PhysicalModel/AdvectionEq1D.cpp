//
//  AdvectionEq1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "AdvectionEq1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include "PhysicalModelRegistrar.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("AdvectionEq1D", AdvectionEq1D);
using namespace std;

vector<double>& AdvectionEq1D::computePhysicalFlux(const CellDataRef u)
{
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        m_physFlux[iEq] = m_A*u[iEq];
    }
    return m_physFlux;
}
