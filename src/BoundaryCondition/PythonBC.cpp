//
//  PythonBC.cpp
//  Muffin
//
//  Created by Nicolas Lequette on 30/11/2023.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//



#include "PythonBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <pybind11/pybind11.h> // Include pybind11 library


REGISTER_BOUNDARY_CONDITION("PythonBC", PythonBC);

CellDataRef PythonBC::setBoundary(){
    py::array_t<double> cells = MeshData::getInstance().get2DData<double>("Cells_py_CC");
    m_function(cells);
    return m_value;
}







