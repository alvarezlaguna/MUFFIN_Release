//
//  FiveMomentQMOM1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 56/09/18.
//  Copyright © 5018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FiveMomentQMOM1D_hpp
#define FiveMomentQMOM1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../CollisionalData/CollisionalData.hpp"
#include "../Parameters.hpp"

using namespace std;
namespace py = pybind11;


class FiveMomentQMOM1D : public PhysicalModel
{
public:
    FiveMomentQMOM1D(string name) : PhysicalModel(name), m_physFlux(NBEQS), m_eigenvals(NBEQS) {
       
    } // The contructor for the case with a single fluid
    ~FiveMomentQMOM1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
};

#endif /* FiveMomentQMOM1D_hpp */
