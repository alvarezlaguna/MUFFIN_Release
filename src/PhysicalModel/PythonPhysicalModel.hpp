//
//  PythonPhysicalModel.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PythonPhysicalModel_hpp
#define PythonPhysicalModel_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../CollisionalData/CollisionalData.hpp"
#include <span>

namespace py = pybind11;

using namespace std;


class PythonPhysicalModel : public PhysicalModel
{
public:
    PythonPhysicalModel(string name) : PhysicalModel(name), m_physFlux(NBEQS), m_eigenvals(NBEQS) {
        m_flux_function     = CollisionalData::getInstance().getData<py::function>("PM_fluxFunction");
        m_maxEigen_function = CollisionalData::getInstance().getData<py::function>("PM_maxEigen_function");

        m_physFlux_view = py::array_t<double>(m_physFlux.size(), m_physFlux.data(),py::cast(m_physFlux));
    }

    ~PythonPhysicalModel(){/*m_physFlux.clear();*/ }
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
protected:
    vector<double> m_eigenvals;
    py::function m_flux_function;
    py::function m_maxEigen_function;
    vector<double> m_physFlux; 
    py::array_t<double> m_physFlux_view;
};

#endif /* PythonPhysicalModel_hpp */
