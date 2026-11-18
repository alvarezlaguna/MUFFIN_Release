//
//  PythonEFieldSourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PythonEFieldSourceTerm_hpp
#define PythonEFieldSourceTerm_hpp

#include <stdio.h>
#include <vector>
#include "SourceTerm.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Parameters;
using namespace std;

class PythonEFieldSourceTerm : public SourceTerm
{
public:
    PythonEFieldSourceTerm(string name) : SourceTerm(name) {
    }
    void setup();
    ~PythonEFieldSourceTerm(){}
    virtual void computeSource();

    private:
        py::function m_function;
};

#endif /* PythonEFieldSourceTerm_hpp */
