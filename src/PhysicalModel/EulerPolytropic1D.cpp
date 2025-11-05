//
//  EulerPolytropic1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 30/09/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "EulerPolytropic1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("EulerPolytropic1D", EulerPolytropic1D);
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& EulerPolytropic1D::getEigenvalues(const CellDataRef u)
{
    const double gamma      = m_polytropicIndex;
    const double R          = m_polytropicConst;
    const double rho        = u[0];
    const double soundSpeed = sqrt(gamma*R*pow(rho, gamma - 1.));
    const double m      = u[1];
    const double v      = m/rho;
    
    m_eigenvals[0] = abs(v - soundSpeed);
    m_eigenvals[1] = abs(v + soundSpeed);
    
    return m_eigenvals;
}

double EulerPolytropic1D::getMaxEigenvalue(const CellDataRef u)
{
    
    m_eigenvals = EulerPolytropic1D::getEigenvalues(u);
    double maxEigenVal = *max_element(m_eigenvals.begin(), m_eigenvals.end());
    double minEigenVal = *min_element(m_eigenvals.begin(), m_eigenvals.end());
    
    return max(abs(maxEigenVal), abs(minEigenVal));
}

vector<double>& EulerPolytropic1D::computePhysicalFlux(const CellDataRef u)
{
    const double rho    = u[0];
    const double m      = u[1];
    m_physFlux[0] = m;
    m_physFlux[1] = m*m/rho + m_polytropicConst*pow(rho,m_polytropicIndex);
        
    return m_physFlux;
}
