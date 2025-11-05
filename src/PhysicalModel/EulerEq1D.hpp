//
//  EulerEq1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef EulerEq1D_hpp
#define EulerEq1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class EulerEq1D : public PhysicalModel
{
public:
    EulerEq1D(string name) : PhysicalModel(name), m_physFlux(3), m_eigenvals(3), m_gamma(GAMMA) {}
    EulerEq1D(string name, double gamma) : PhysicalModel(name), m_physFlux(3), m_eigenvals(3), m_gamma(gamma) {}

    ~EulerEq1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    double getGamma(){return m_gamma;}
protected:
    vector<double> m_physFlux;    
    vector<double> m_eigenvals;
    double m_gamma;
};

#endif /* EulerEq1D_hpp */
