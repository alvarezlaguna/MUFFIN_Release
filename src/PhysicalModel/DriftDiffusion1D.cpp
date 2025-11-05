//
//  DriftDiffusion1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DriftDiffusion1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>
using namespace std;

REGISTER_PHYSICALMODEL("DriftDiffusion1D", DriftDiffusion1D);

vector<double>& DriftDiffusion1D::computePhysicalFlux(const CellDataRef u)
{
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        m_physFlux[iEq] = 0.;
    }
    return m_physFlux;
}
