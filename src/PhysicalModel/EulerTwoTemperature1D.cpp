//
//  EulerTwoTemperature1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 03/10/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "EulerTwoTemperature1D.hpp"
#include <vector>
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

#include "PhysicalModelRegistrar.hpp"

REGISTER_PHYSICALMODEL("EulerTwoTemperature1D", EulerTwoTemperature1D);


using namespace std;


vector<double>& EulerTwoTemperature1D::getEigenvalues(const CellDataRef u)
{
    const double gamma  = 3.;
    const double rho    = u[0];
    const double m      = u[1];
    const double ex     = u[2];
    const double v      = m/rho;
    const double px     = ex - m*m/(rho);
    double soundSpeed = sqrt(gamma*(ex - m*m/(rho))/rho);
    
    m_eigenvals[0] = (v - soundSpeed);
    m_eigenvals[1] = (v);
    m_eigenvals[2] = (v);
    m_eigenvals[3] = (soundSpeed + v );
    
    return m_eigenvals;
}

double EulerTwoTemperature1D::getMaxEigenvalue(const CellDataRef u)
{
    
    m_eigenvals = EulerTwoTemperature1D::getEigenvalues(u);
    double maxEigenVal = *max_element(m_eigenvals.begin(), m_eigenvals.end());
    double minEigenVal = *min_element(m_eigenvals.begin(), m_eigenvals.end());
    
    //cout<<"maxEigenVal = "<<maxEigenVal<<"\n";
    //cout<<"minEigenVal = "<<minEigenVal<<"\n";
    
    return max(abs(maxEigenVal), abs(minEigenVal));
}

vector<double>& EulerTwoTemperature1D::computePhysicalFlux(const CellDataRef u)
{
    double rho          = u[0];
    const double m      = u[1];
    const double ex     = u[2];
    const double py     = u[3];

    const double v      = m/rho;
    const double px     = ex - m*m/(rho);

    m_physFlux[0] = m;
    m_physFlux[1] = ex;
    m_physFlux[2] = rho*v*v*v + 3*px*v;
    m_physFlux[3] = py*v;
    
    return m_physFlux;
}
