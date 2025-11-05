//
//  EulerEq1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "EulerEq1D.hpp"
#include "PhysicalModelRegistrar.hpp"
#include "PhysicalModelRegistrar.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("EulerEq1D", EulerEq1D);
#include <cmath>
using namespace std;


vector<double>& EulerEq1D::getEigenvalues(const CellDataRef u)
{
    const double gamma  = GAMMA;
    const double rho    = u[0];
    const double m      = u[1];
    const double e      = u[2];
    const double v      = m/rho;
    const double soundSpeed = sqrt(gamma*(gamma - 1)*(e - m*m/(2*rho))/rho);

    m_eigenvals[0] = v - soundSpeed;
    m_eigenvals[1] = v;
    m_eigenvals[2] = v + soundSpeed;
    
    return m_eigenvals;
}

double EulerEq1D::getMaxEigenvalue(const CellDataRef u)
{
    const double gamma  = GAMMA;
    const double rho    = u[0];
    const double m      = u[1];
    const double e      = u[2];
    const double v      = m/rho;
    const double soundSpeed = sqrt(gamma*(gamma - 1)*(e - m*m/(2*rho))/rho);
    
    m_eigenvals[0] = v - soundSpeed;
    m_eigenvals[1] = v;
    m_eigenvals[2] = v + soundSpeed;
    
    return max(abs(m_eigenvals[0]), max(abs(m_eigenvals[1]), abs(m_eigenvals[2])));
}

vector<double>& EulerEq1D::computePhysicalFlux(const CellDataRef u)
{
    const double gamma  = GAMMA;
    const double rho    = u[0];
    const double m      = u[1];
    const double e      = u[2];
    m_physFlux[0] = m;
    m_physFlux[1] = m*m/rho*(3 - gamma)/2 + (gamma - 1)*e;
    m_physFlux[2] = e*m/rho*gamma - (gamma - 1)/2*pow(m,3)/pow(rho,2);

    return m_physFlux;
}
