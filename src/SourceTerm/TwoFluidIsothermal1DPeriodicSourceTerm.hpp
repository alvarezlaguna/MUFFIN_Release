//
//  TwoFluidIsothermal1DPeriodicSourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 06/02/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermal1DPeriodicSourceTerm_hpp
#define TwoFluidIsothermal1DPeriodicSourceTerm_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "SourceTerm.hpp"
#include "../LinearSolver/LinearSolver.hpp"

using namespace Parameters;
using namespace std;

typedef vector< vector<double> > Matrix;

class TwoFluidIsothermal1DPeriodicSourceTerm : public SourceTerm
{
public:
    TwoFluidIsothermal1DPeriodicSourceTerm(string name) : SourceTerm(name), m_A(NBCELLS,vector<double>(NBCELLS,0.)), m_B(NBCELLS, 0.)
    {
        // Set up initial field of Phi
        m_x = Parameters::PHIINITIAL;
        unsigned int n1 = NBCELLS - 1;
        
        // Set up the matrix
        m_A[0][0] = -2;
        m_A[0][1] = 1;
        m_A[0][n1] = 1;
        // Inner cells
        for (unsigned int i = 1; i < n1; i++) {
            m_A[i][i]     = -2;
            m_A[i][i - 1] = 1;
            m_A[i][i + 1] = 1;
        }
        // Last cell
        m_A[n1][0] = 1;
        m_A[n1][n1]   = -2;
        m_A[n1][n1-1] = 1;
    }
    ~TwoFluidIsothermal1DPeriodicSourceTerm(){}
    void setup();
    virtual void computeSource();
    double computeIonizationConstant();
    
private:
    /// TODO: write here the matrix and the vectors for the linear system
    unique_ptr<LinearSolver> m_linearSolver;
    Matrix m_A;         // Matrix of the linear solver
    vector<double> m_x; // Solution of the linear solver
    vector<double> m_B; // rhs of the linear solver
    
};


#endif /* TwoFluidIsothermal1DPeriodicSourceTerm_hpp */
