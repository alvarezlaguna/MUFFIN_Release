//
//  FiveMomentQMOM1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FiveMomentQMOM1D.hpp"
#include <vector>
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

#include "PhysicalModelRegistrar.hpp"

REGISTER_PHYSICALMODEL("FiveMomentQMOM1D", FiveMomentQMOM1D);

using namespace std;


vector<double>& FiveMomentQMOM1D::getEigenvalues(const CellDataRef u)
{
    // We compute just the max
    m_eigenvals[0] = FiveMomentQMOM1D::getMaxEigenvalue(u);
    return m_eigenvals;
}

double FiveMomentQMOM1D::getMaxEigenvalue(const CellDataRef u)
{

    const double rho    = u[0];
    const double v      = u[1]/rho;
    const double p      = u[2] - rho*v*v;
    const double q      = u[3] - rho*v*v*v - 3*v*p;
    const double r      = u[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;
    const double r_Max  =  3*p*p/rho;
    const double K      = r - r_Max;
    const double a_n     = sqrt(p/rho);
    const double q_nStar = q/(rho*pow(a_n,3));
    const double K_nStar = K/(rho*pow(a_n,4));
    const double r_nStar = K_nStar + 3.;

    const double xi_nStar = max(r_nStar - pow(q_nStar,2) - 1., 0.); // It should be positive, otherwise it is numerical error
    if(r_nStar - pow(q_nStar,2) - 1. < 0.){
        py::print("WARNING: You are out of the realizability");
        py::print("QStar = ", q_nStar, "  RStar = ", r_nStar);
    }
    const double A = q_nStar/2.;
    const double B = 1. + xi_nStar + pow(q_nStar,2)/4.;
    const double C = sqrt(xi_nStar*(1 + xi_nStar));

    return max(abs(v + a_n*(A + sqrt(B + C))), max(abs(v + a_n*(A + sqrt(B - C))), max(abs(v + a_n*(A - sqrt(B + C))), abs(v + a_n*(A - sqrt(B - C))))));
}

vector<double>& FiveMomentQMOM1D::computePhysicalFlux(const CellDataRef u)
{
    const double rho    = u[0];
    const double v      = u[1]/rho;
    const double p      = u[2] - rho*v*v;
    const double q      = u[3] - rho*v*v*v - 3*v*p;
    const double r      = u[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;
    const double r_Max  =  3*p*p/rho;
    const double K      = r - r_Max;
    m_physFlux[0] = u[1];
    m_physFlux[1] = u[2];
    m_physFlux[2] = u[3];
    m_physFlux[3] = u[4];
    m_physFlux[4] = rho*v*v*v*v*v + 10.*v*v*v*p + 10.*v*v*q + 5.*v*r;
    
    m_physFlux[4] += 2*r*q/p - pow(q,3)/pow(p,2);
    // cout<<"s = "<< s.cast<double>() <<"\n";
    
    return m_physFlux;
}
