//
//  PythonBC.hpp
//  Muffin
//
//  Created by Nicolas Lequette on 30/11/2023.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PythonBC_hpp
#define PythonBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../CollisionalData/CollisionalData.hpp"
namespace py = pybind11;

class PythonBC : public BoundaryCondition
{
public:
    PythonBC(string name, string side) : BoundaryCondition(name, side){ 
        if (side == "Left"){
            m_function = CollisionalData::getInstance().getData<py::function>("functionInlet");

        }
        if (side == "Right"){
            m_function = CollisionalData::getInstance().getData<py::function>("functionOutlet");
        }
        
      }
    ~PythonBC() {}
    virtual CellDataRef setBoundary();
private:
    py::function m_function;
};

#endif /* PythonBC_hpp */
