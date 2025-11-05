//
//  AdvectionEq1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef AdvectionEq1D_hpp
#define AdvectionEq1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class AdvectionEq1D : public PhysicalModel
{
public:
    AdvectionEq1D(string name) : PhysicalModel(name), m_physFlux(Parameters::NBEQS), m_eigenvals(Parameters::NBEQS) {}
    ~AdvectionEq1D(){m_physFlux.clear();}
    virtual vector<double>& getEigenvalues(const CellDataRef u)
        { for (int iEq = 0; iEq < Parameters::NBEQS; iEq++){m_eigenvals[iEq] = m_A;} return m_eigenvals;}
    virtual double getMaxEigenvalue(const CellDataRef u){return m_A;}
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
protected:
    const double m_A = Parameters::A;
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    
};

#endif /* AdvectionEq1D_hpp */
