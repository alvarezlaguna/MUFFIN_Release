//
//  TimeMethod.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TimeMethod_hpp
#define TimeMethod_hpp

#include <stdio.h>
#include "../Component.hpp"
#include "../SpaceMethod/SpaceMethod.hpp"
#include "../DiffusiveFlux/DiffusiveFlux.hpp"
#include <iostream>

class TimeMethod : public Component
{
public:
    TimeMethod(string name) : Component(name), m_iter(0), m_norm(NBEQS, 1e10){}
    virtual ~TimeMethod(){}
    virtual void setup() {
        m_cells  = &MeshData::getInstance().getData<Cell1D>("Cells");
        m_rhs    = &MeshData::getInstance().getData<double>("rhs");
        m_source = &MeshData::getInstance().get2DData<double>("source");

    }

    virtual void takeStep(double dt) = 0;
    void takeStep(){
        m_spaceMethod->discretize();
        double dt = m_spaceMethod->getDt();
        this->takeStep(dt);
        }

    void setSpaceMethod(SpaceMethod* sp){m_spaceMethod = sp;}
    void setDiffusiveFlux(DiffusiveFlux* df){m_diffFlux = df;}
    vector<double>& getNorm(){return m_norm;}
    unsigned int getIter(){return m_iter;}
    
protected:
    unsigned int m_iter;
    vector<double> m_norm;
    SpaceMethod* m_spaceMethod;
    DiffusiveFlux* m_diffFlux;
    
    vector<Cell1D>* m_cells; // Pointer to the cells to avoid lookup
    vector<double>* m_rhs; // Pointer to the cells to avoid lookup
    py::array_t<double>* m_source; // Pointer to the cells to avoid lookup

};

#endif /* TimeMethod_hpp */
