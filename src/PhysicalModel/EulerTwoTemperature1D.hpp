//
//  EulerTwoTemperature1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 03/10/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef EulerTwoTemperature1D_hpp
#define EulerTwoTemperature1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class EulerTwoTemperature1D : public PhysicalModel
{
public:
    EulerTwoTemperature1D(string name) : PhysicalModel(name), m_physFlux(4), m_eigenvals(4) {} // The contructor for the case with a single fluid
    ~EulerTwoTemperature1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    //double getSoundSpeed(){return m_soundSpeed;}
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;

};

#endif /* EulerTwoTemperature1D_hpp */