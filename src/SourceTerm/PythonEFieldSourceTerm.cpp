//
//  PythonEFieldSourceTerm.cpp
//  Muffin
//
//  Created by Nicolas Lequette and Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include <chrono>
#include "PythonEFieldSourceTerm.hpp"

// Self-register with the factory
#include "SourceTermRegistrar.hpp"
REGISTER_SOURCETERM("PythonEFieldSourceTerm", PythonEFieldSourceTerm);
#include "../MeshData/MeshData.hpp"
#include "../CollisionalData/CollisionalData.hpp"
#include <stdlib.h>     //for using the function sleep


namespace py = pybind11;

using namespace std;

void PythonEFieldSourceTerm::setup()
{

    SourceTerm::setup();


    m_function  = CollisionalData::getInstance().getData<py::function>("function_EFieldPhi");
}

void PythonEFieldSourceTerm::computeSource()
{
    py::array_t<double> EFieldPhi = MeshData::getInstance().get2DData<double>("EFieldPhi");
    // old boundary call  = MeshData::getInstance().getData<double>("boundaries");


    // py::capsule free_when_done(m_ptr, [](void *f) {
    //    double *foo = reinterpret_cast<double *>(f);
    //    delete[] foo;
    // });

    py::array_t<double> m_pythonArray=MeshData::getInstance().get2DData<double>("Cells_py_CC");
           

    py::object nu=m_function(m_pythonArray, EFieldPhi);
    


    /* Getting number of milliseconds as a double. */







    // py::print(m_pythonArray);
    // sleep(2);
    // exit(0);
}
