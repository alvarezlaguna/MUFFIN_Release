//
//  EulerPolytropic1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 30/09/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef EulerPolytropic1D_hpp
#define EulerPolytropic1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"

using namespace std;


class EulerPolytropic1D : public PhysicalModel
{
public:
    EulerPolytropic1D(string name) : PhysicalModel(name), m_physFlux(2), m_eigenvals(2), m_polytropicIndex(POLYTROPICINDEX[0]), m_polytropicConst(1.) {} // The contructor for the case with a single fluid
    EulerPolytropic1D(string name, double polytropicIndex, double polytropicConst) : PhysicalModel(name), m_physFlux(2), m_eigenvals(2), m_polytropicIndex(polytropicIndex), m_polytropicConst(polytropicConst){}
    ~EulerPolytropic1D(){m_physFlux.clear();}
    
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    //double getSoundSpeed(){return m_soundSpeed;}
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    double m_polytropicIndex;
    double m_polytropicConst;
};

#endif /* EulerPolytropic1D_hpp */
