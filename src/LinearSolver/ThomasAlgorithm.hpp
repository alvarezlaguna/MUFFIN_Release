//
//  ThomasAlgorithm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 24/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ThomasAlgorithm_hpp
#define ThomasAlgorithm_hpp

#include <stdio.h>
#include "LinearSolver.hpp"

class ThomasAlgorithm : public LinearSolver
{
public:
    ThomasAlgorithm(string name, unsigned int size) : LinearSolver(name, size), m_cprime(size, 0.), m_dprime(size, 0){}
    ~ThomasAlgorithm(){}
    void setup(){}
    virtual void solveLinearSystem(const Matrix& A, vector<double>& x, const vector<double>& B);
protected:
    vector<double> m_cprime;
    vector<double> m_dprime;
};


#endif /* ThomasAlgorithm_hpp */
