//
//  LinearizedGrad1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LinearizedGrad1D_hpp
#define LinearizedGrad1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class LinearizedGrad1D : public PhysicalModel
{
public:
    LinearizedGrad1D(string name) : PhysicalModel(name), m_physFlux(Parameters::LEVEL), m_eigenvals(Parameters::LEVEL) {
        m_eigenvals = Parameters::EIGENVALUES_GRAD;
        double maxEigenVal = *max_element(m_eigenvals.begin(), m_eigenvals.end());
        double minEigenVal = *min_element(m_eigenvals.begin(), m_eigenvals.end());
        m_maxEigenvalue = max(abs(maxEigenVal), abs(minEigenVal));
    } // The contructor for the case with a single fluid
    ~LinearizedGrad1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);

    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    double m_maxEigenvalue;
};

#endif /* LinearizedGrad1D_hpp */
