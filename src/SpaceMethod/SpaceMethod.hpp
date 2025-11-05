//
//  SpaceMethod.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 16/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SpaceMethod_hpp
#define SpaceMethod_hpp

#include <stdio.h>
#include <iostream>
#include "../Component.hpp"
#include "../PhysicalModel/PhysicalModel.hpp"
#include <mpi.h>

using namespace std;

class SpaceMethod : public Component
{
public:
    SpaceMethod(string name) : Component(name), m_dtOvdx(0) {}
    void discretize(){
        setBoundaries();
        computeRHS();
    }
    virtual ~SpaceMethod(){}
    void setPhysicalModel(PhysicalModel* pm){m_pm = pm;}
    PhysicalModel* getPhysicalModel(){return m_pm;}
    virtual double getDtOvDx() = 0;
    virtual double getDt() = 0;
    virtual void setBoundaries() = 0;
    virtual void computeRHS() = 0;
protected:
    PhysicalModel* m_pm;
    double m_dtOvdx;
    double m_dt;
};

#endif /* SpaceMethod_hpp */
