//
//  ThomasAlgorithmPeriodic.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 30/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ThomasAlgorithmPeriodic.hpp"
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void ThomasAlgorithmPeriodic::solveLinearSystem(const Matrix& A, vector<double>& x, const vector<double>& B)
{
    py::module numpy = py::module::import("numpy");
    py::module linalg = py::module::import("numpy.linalg");
    
    py::buffer_info buf_a = m_a.request();
    py::buffer_info buf_b = m_b.request();
    double *ptra_py = (double *) buf_a.ptr;
    double *ptrb_py = (double *) buf_b.ptr;

    size_t X = buf_b.shape[0];
    size_t Y = buf_a.shape[1];

    // We solve the linear system Ax = b of a cyclic matrix
    // We impose as condition, the average charge to be zero in the domain
    for (size_t idx = 0; idx < X; idx++){
        if (idx != NBCELLS){
            ptrb_py[idx] = B[idx];
            for (size_t idy = 0; idy < Y; idy++){
                if (idy != NBCELLS){
                    ptra_py[idx*Y + idy] = A[idx][idy];
                }
                else{
                    ptra_py[idx*Y + idy] = 1;
                }
            }

        }
        else{
            ptrb_py[idx] = 0.;
            for (size_t idy = 0; idy < Y; idy++){
                if (idy != NBCELLS){
                    ptra_py[idx*Y + idy] = 1.;
                }
            }
        }
    }

    m_x_py = linalg.attr("solve")(m_a, m_b);
    py::buffer_info buf_x = m_x_py.request();
    double *ptrx_py = (double *) buf_x.ptr;
    
    for (size_t idx = 0; idx < NBCELLS; idx++){
        x[idx] = ptrx_py[idx];
    }
        
//    cout<<"A = \n";
//    for(unsigned int i = 0; i < NBCELLS; i++){
//        cout<<"| ";
//        for(unsigned int j = 0; j < NBCELLS; j++){
//            cout<< A[i][j] <<"\t";
//        }
//        cout<<"| \n";
//    }
//    cout<<"B = \n";
//    for(unsigned int i = 0; i < NBCELLS; i++){
//        cout<<"| ";
//        cout<< B[i] <<"\t";
//        cout<<"| \n";
//    }
    
//    // First system
//    m_cprime[0] = A[0][1]/(2*A[0][0]);
//    m_dprime[0] = B[0]/(2*A[0][0]);
//    // Forward sweep
//    unsigned int n1 = NBCELLS - 1;
//    for(unsigned int i = 1; i < n1; i++)
//    {
//        double a_i        = A[i][i-1];
//        double b_i        = A[i][i];
//        double c_i        = A[i][i+1];
//        double d_i        = B[i];
//        m_cprime[i]  = c_i/(b_i - a_i*m_cprime[i - 1]);
//        m_dprime[i]  = (d_i - a_i*m_dprime[i - 1])/(b_i - a_i*m_cprime[i - 1]);
//    }
//    // Last cell
//    double a_i        = A[n1][n1-1];
//    double b_i        = A[n1][n1] + A[n1][0]*A[0][n1]/A[0][0];
//    double d_i        = B[n1];
//    m_dprime[n1]  = (d_i - a_i*m_dprime[n1 - 1])/(b_i - a_i*m_cprime[n1 - 1]);
//    // Solution
//    m_y[n1] = m_dprime[n1];
//    for (int i = n1 - 1; i >= 0; i--) {
//        m_y[i] = m_dprime[i] - m_cprime[i]*m_y[i + 1];
//    }
//    
////    cout<<"y = \n";
////    for(unsigned int i = 0; i < NBCELLS; i++){
////        cout<<"| ";
////        cout<< m_y[i] <<"\t";
////        cout<<"| \n";
////    }
//    
//    // second system
//    m_cprime[0] = A[0][1]/(2*A[0][0]);
//    m_dprime[0] = -A[0][0]/(2*A[0][0]);
//    // Forward sweep
//
//    for(unsigned int i = 1; i < n1; i++)
//    {
//        double a_i        = A[i][i-1];
//        double b_i        = A[i][i];
//        double c_i        = A[i][i+1];
//        double d_i        = 0.;
//        m_cprime[i]  = c_i/(b_i - a_i*m_cprime[i - 1]);
//        m_dprime[i]  = (d_i - a_i*m_dprime[i - 1])/(b_i - a_i*m_cprime[i - 1]);
//    }
//    // Last cell
//    a_i        = A[n1][n1-1];
//    b_i        = A[n1][n1] + A[n1][0]*A[0][n1]/A[0][0];
//    d_i        = A[n1][0];
//    m_dprime[n1]  = (d_i - a_i*m_dprime[n1 - 1])/(b_i - a_i*m_cprime[n1 - 1]);
//    
//    // Solution
//    m_q[n1] = m_dprime[n1];
//    for (int i = n1 - 1; i >= 0; i--) {
//        m_q[i] = m_dprime[i] - m_cprime[i]*m_q[i + 1];
//    }
////    cout<<"q = \n";
////    for(unsigned int i = 0; i < NBCELLS; i++){
////        cout<<"| ";
////        cout<< m_q[i] <<"\t";
////        cout<<"| \n";
////    }
//    
//    double a1ovb1 = A[0][n1]/A[0][0];
//    double factor = (m_y[0] - a1ovb1*m_y[n1])/(1 + (m_q[0] - a1ovb1*m_q[n1]));
//    
//    //cout<<"factor = "<< factor <<"\n";
//    
//    //cout<<"\nx = \n";
//    for(unsigned int i = 0; i < m_y.size(); i++)
//    {
//        x[i] = m_y[i] - factor*m_q[i];
//    }
//    //cout<<"exit here \n";
//    //exit(0);
}
