//
//  EulerIsothermal1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "EulerIsothermal1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("EulerIsothermal1D", EulerIsothermal1D);
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& EulerIsothermal1D::getEigenvalues(const CellDataRef u)
{
    const double soundSpeed = m_soundSpeed;
    const double rho    = u[0];
    const double m      = u[1];
    const double v      = m/rho;
    
    m_eigenvals[0] = v - soundSpeed;
    m_eigenvals[1] = v + soundSpeed;
    
    return m_eigenvals;
}

double EulerIsothermal1D::getMaxEigenvalue(const CellDataRef u)
{

    m_eigenvals = EulerIsothermal1D::getEigenvalues(u);
    double maxEigenVal = *max_element(m_eigenvals.begin(), m_eigenvals.end());
    double minEigenVal = *min_element(m_eigenvals.begin(), m_eigenvals.end());
    
    return max(abs(maxEigenVal), abs(minEigenVal));
}

vector<double>& EulerIsothermal1D::computePhysicalFlux(const CellDataRef u)
{
    const double soundSpeed = m_soundSpeed;
    const double rho    = u[0];
    const double m      = u[1];
    m_physFlux[0] = m;
    m_physFlux[1] = m*m/rho + rho*soundSpeed*soundSpeed;
    
    return m_physFlux;
}
