//
//  PythonSourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PythonSourceTerm_hpp
#define PythonSourceTerm_hpp

#include <stdio.h>
#include <vector>
#include "SourceTerm.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Parameters;
using namespace std;

class PythonSourceTerm : public SourceTerm
{
public:
    PythonSourceTerm(string name) : SourceTerm(name) {
    }
    void setup();
    ~PythonSourceTerm(){}
    virtual void computeSource();

    private:
        py::function m_function;
};

#endif /* PythonSourceTerm_hpp */
