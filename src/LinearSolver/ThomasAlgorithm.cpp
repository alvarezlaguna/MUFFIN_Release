//
//  ThomasAlgorithm.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 24/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ThomasAlgorithm.hpp"
#include <iostream>

void ThomasAlgorithm::solveLinearSystem(const Matrix& A, vector<double>& x, const vector<double>& B)
{
    m_cprime[0] = A[0][1]/A[0][0];
    m_dprime[0] = B[0]/A[0][0];
    // Forward sweep
    unsigned int n1 = NBCELLS - 1;
    for(unsigned int i = 1; i < n1; i++)
    {
        double a_i        = A[i][i-1];
        double b_i        = A[i][i];
        double c_i        = A[i][i+1];
        double d_i        = B[i];
        m_cprime[i]  = c_i/(b_i - a_i*m_cprime[i - 1]);
        m_dprime[i]  = (d_i - a_i*m_dprime[i - 1])/(b_i - a_i*m_cprime[i - 1]);
    }
    // Last cell
    double a_i        = A[n1][n1-1];
    double b_i        = A[n1][n1];
    double d_i        = B[n1];
    m_dprime[n1]  = (d_i - a_i*m_dprime[n1 - 1])/(b_i - a_i*m_cprime[n1 - 1]);
    // Solution
    x[n1] = m_dprime[n1];
    for (int i = n1 - 1; i >= 0; i--) {
        x[i] = m_dprime[i] - m_cprime[i]*x[i + 1];
    }
}
