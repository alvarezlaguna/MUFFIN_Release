//
//  SourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SourceTerm_hpp
#define SourceTerm_hpp

#include <stdio.h>
#include "../Component.hpp"
#include "../PhysicalModel/PhysicalModel.hpp"


class SourceTerm : public Component
{
public:
    SourceTerm(string name) : Component(name), m_frequency(0.) {}
    virtual ~SourceTerm(){}
    virtual void setup(){
        m_cells = &MeshData::getInstance().getData<Cell1D>("Cells");
    };
    void setPhysicalModel(PhysicalModel* pm){m_pm = pm;}
    virtual void computeSource() = 0;
    void setFrequency(double freq){m_frequency = freq;};
    double getFrequency(){return m_frequency;}
    
protected:
    PhysicalModel* m_pm;
    vector<Cell1D>* m_cells; // Pointer to the cells to avoid lookup

private:
    double m_frequency;
};

#endif /* SourceTerm_hpp */
