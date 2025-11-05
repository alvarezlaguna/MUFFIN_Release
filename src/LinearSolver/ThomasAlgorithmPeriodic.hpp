//
//  ThomasAlgorithmPeriodic.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 30/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ThomasAlgorithmPeriodic_hpp
#define ThomasAlgorithmPeriodic_hpp

#include <stdio.h>
#include "LinearSolver.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class ThomasAlgorithmPeriodic : public LinearSolver
{
public:
    ThomasAlgorithmPeriodic(string name, unsigned int size) : LinearSolver(name, size), m_b(NBCELLS + 1), m_x_py(NBCELLS + 1){
        //Initialize A
        py::module numpy = py::module::import("numpy");
        py::module linalg = py::module::import("numpy.linalg");
        
        py::tuple args = py::make_tuple(NBCELLS + 1, NBCELLS + 1);
        m_a = numpy.attr("zeros")(args);

    }
    ~ThomasAlgorithmPeriodic(){}
    void setup(){}
    virtual void solveLinearSystem(const Matrix& A, vector<double>& x, const vector<double>& B);
protected:
    py::array_t<double> m_b;
    py::array_t<double> m_x_py;
    py::array_t<double> m_a;
};


#endif /* ThomasAlgorithmPeriodic_hpp */
