//
//  DriftDiffusion1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DriftDiffusion1D_hpp
#define DriftDiffusion1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class DriftDiffusion1D : public PhysicalModel
{
public:
    DriftDiffusion1D(string name) : PhysicalModel(name), m_physFlux(Parameters::NBEQS), m_eigenvals(Parameters::NBEQS) {}
    ~DriftDiffusion1D(){m_physFlux.clear();}
    virtual vector<double>& getEigenvalues(const CellDataRef u)
        { for (int iEq = 0; iEq < Parameters::NBEQS; iEq++){m_eigenvals[iEq] = 1.;} return m_eigenvals;}
    virtual double getMaxEigenvalue(const CellDataRef u){return 1.;}
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    
};

#endif /* DriftDiffusion1D_hpp */
