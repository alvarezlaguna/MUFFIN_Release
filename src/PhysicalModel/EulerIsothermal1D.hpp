//
//  EulerIsothermal1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef EulerIsothermal1D_hpp
#define EulerIsothermal1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class EulerIsothermal1D : public PhysicalModel
{
public:
    EulerIsothermal1D(string name) : PhysicalModel(name), m_physFlux(2), m_eigenvals(2), m_soundSpeed(SOUNDSPEED[0]) {} // The contructor for the case with a single fluid
    EulerIsothermal1D(string name, double soundSpeed) : PhysicalModel(name), m_physFlux(2), m_eigenvals(2), m_soundSpeed(soundSpeed){}
    ~EulerIsothermal1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    double getSoundSpeed(){return m_soundSpeed;}
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    double m_soundSpeed;
};

#endif /* EulerIsothermal1D_hpp */
