//
//  LinearizedGrad1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "LinearizedGrad1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("LinearizedGrad1D", LinearizedGrad1D);
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& LinearizedGrad1D::getEigenvalues(const CellDataRef u)
{
    return m_eigenvals;    
}

double LinearizedGrad1D::getMaxEigenvalue(const CellDataRef u)
{
    return m_maxEigenvalue;
}

vector<double>& LinearizedGrad1D::computePhysicalFlux(const CellDataRef u)
{
    m_physFlux[0] = u[1];
    m_physFlux[Parameters::LEVEL - 1] = sqrt(Parameters::LEVEL) * u[Parameters::LEVEL - 2];
    for (int i = 1; i < Parameters::LEVEL - 1; ++i){
        m_physFlux[i] = sqrt(i+1) * u[i+1] + sqrt(i) * u[i-1];
    }
    return m_physFlux;

}
