//
//  TwoFluidIsothermal1DSourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermal1DSourceTerm_hpp
#define TwoFluidIsothermal1DSourceTerm_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "SourceTerm.hpp"
#include "../LinearSolver/LinearSolver.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../MeshData/MeshData.hpp"

using namespace Parameters;
using namespace std;

typedef vector< vector<double> > Matrix;

class TwoFluidIsothermal1DSourceTerm : public SourceTerm
{
public:
    TwoFluidIsothermal1DSourceTerm(string name) : SourceTerm(name), m_A(NBCELLS,vector<double>(NBCELLS,0.)), m_B(NBCELLS, 0.)
    {
        // Set up initial field of Phi
        m_x = Parameters::PHIINITIAL;
        vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
        vector<double>& x = MeshData::getInstance().getData<double>("x");
        // Implementation with non-homogeneous mesh
        double x_imO = x[0] - cells[0].dx;
        double x_i   = x[0];
        double x_ipO = x[1];
        double denominator  = (x_ipO - x_imO)*(x_ipO - x_i)*(x_i - x_imO);
        
        
        m_A[0][0]  = -2*(x_ipO - x_imO)/denominator;
        m_A[0][1] = 2*(x_i - x_imO)/denominator;
        // Inner cells
        unsigned int n1 = NBCELLS - 1;
        for (unsigned int i = 1; i < n1; i++) {
            x_imO = x[i - 1];
            x_i   = x[i];
            x_ipO = x[i + 1];
            denominator  = (x_ipO - x_imO)*(x_ipO - x_i)*(x_i - x_imO);
            m_A[i][i]     = -2*(x_ipO - x_imO)/denominator;
            m_A[i][i - 1] = 2*(x_ipO - x_i)/denominator;
            m_A[i][i + 1] = 2*(x_i - x_imO)/denominator;
        }
        // Last cell
        x_imO = x[n1 - 1];
        x_i   = x[n1];
        x_ipO = x[n1] + cells[n1].dx;
        m_A[n1][n1]   = -2*(x_ipO - x_imO)/denominator;
        m_A[n1][n1-1] = 2*(x_ipO - x_i)/denominator;
        
//        cout<<"A[0][0] = "<<m_A[0][0]<<"\t A[0][1] = "<<m_A[0][1]<<"\n";
//        for(int iCell = 1; iCell < n1; iCell++){
//            cout<<"A["<<iCell<<"]["<<iCell - 1<<"] = "<<m_A[iCell][iCell-1]<<"\t A["<<iCell<<"]["<<iCell <<"] = "<<m_A[iCell][iCell]<<"\tA["<<iCell<<"]["<<iCell + 1<<"] = "<<m_A[iCell][iCell+1]<<"\n";
//        }
//        cout<<"A["<<n1<<"]["<<n1 -1<<"] = "<<m_A[n1][n1 - 1]<<"\t A["<<n1<<"]["<<n1<<"] = "<<m_A[n1][n1]<<"\n";
//        exit(0);

        
        
//        // Set up the matrix
//        m_A[0][0] = -2;
//        m_A[0][1] = 1;
//        unsigned int n1 = NBCELLS - 1;
//        // Inner cells
//        for (unsigned int i = 1; i < n1; i++) {
//            m_A[i][i]     = -2;
//            m_A[i][i - 1] = 1;
//            m_A[i][i + 1] = 1;
//        }
//        // Last cell
//        m_A[n1][n1]   = -2;
//        m_A[n1][n1-1] = 1;
    }
    ~TwoFluidIsothermal1DSourceTerm(){}
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

#endif /* TwoFluidIsothermal1DSourceTerm_hpp */
