//
//  PythonPhysicalModel.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "PythonPhysicalModel.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("PythonPhysicalModel", PythonPhysicalModel);
#include <cmath>
using namespace std;

namespace py = pybind11;

vector<double> &PythonPhysicalModel::getEigenvalues(const CellDataRef u)
{
    const double gamma = GAMMA;
    const double rho = u[0];
    const double m = u[1];
    const double e = u[2];
    const double v = m / rho;
    const double soundSpeed = sqrt(gamma * (gamma - 1) * (e - m * m / (2 * rho)) / rho);

    m_eigenvals[0] = v - soundSpeed;
    m_eigenvals[1] = v;
    m_eigenvals[2] = v + soundSpeed;

    return m_eigenvals;
}

double PythonPhysicalModel::getMaxEigenvalue(const CellDataRef u)
{

    py::object max_eigen = m_maxEigen_function(u.getView());

    return max_eigen.cast<double>();
    
}

vector<double> &PythonPhysicalModel::computePhysicalFlux(const CellDataRef u)
{


    m_flux_function(u.getView(), m_physFlux_view);

    // cout<<"Cell = "<<u.getCellId()<<"\n";
    // cout<<"size of flux"<<m_physFlux.size()<<"\n";
    // cout<<"u = "<<(u.getCellVector())[0]<<",\t"<< (u.getCellVector())[1] <<",\t"<< (u.getCellVector())[2] <<"\n";
    // cout<<"pysical_flux = "<<m_physFlux[0]<<",\t"<< m_physFlux[1] <<",\t"<< m_physFlux[2] <<"\n";
    // for (unsigned int iEq = 0 ; iEq < m_physFlux.size(); iEq++){
    //     m_physFlux[iEq] = m_physFlux_view.at(iEq);
    // }
    return m_physFlux;
}
