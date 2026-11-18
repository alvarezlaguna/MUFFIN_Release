//
//  SixMomentAQMOM1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 56/09/18.
//  Copyright © 5018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SixMomentAQMOM1D_hpp
#define SixMomentAQMOM1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../CollisionalData/CollisionalData.hpp"
#include "../Parameters.hpp" 

using namespace std;
namespace py = pybind11;


class SixMomentAQMOM1D : public PhysicalModel
{
public:
    SixMomentAQMOM1D(string name) : PhysicalModel(name), m_physFlux(6), m_eigenvals(6) {
       
    } // The contructor for the case with a single fluid
    ~SixMomentAQMOM1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
};

#endif /* SixMomentAQMOM1D_hpp */
