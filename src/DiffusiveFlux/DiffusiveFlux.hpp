//
//  DiffusiveFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DiffusiveFlux_hpp
#define DiffusiveFlux_hpp

#include <stdio.h>
#include "../Component.hpp"
#include "../MeshData/MeshData.hpp"
#include "../CollisionalData/Mixture.hpp"
#include "../MeshData/Cell1D.hpp"


using namespace Parameters;
using namespace std;

class DiffusiveFlux : public Component
{
public :
    DiffusiveFlux(string name) : Component(name) {}
    virtual ~DiffusiveFlux(){}
    
    void setMixture(Mixture* mixture) {m_mixture = mixture;}
    Mixture* getMixture(){return m_mixture;}
    virtual void setup() {
        m_cells = &MeshData::getInstance().getData<Cell1D>("Cells");
    }
    virtual void unsetup() {}
    
    virtual void computeDiffusiveStep(const double dt) {}

protected :
    vector<Cell1D>* m_cells; // Pointer to the cells to avoid lookup
    Mixture* m_mixture;
};
#endif /* DiffusiveFlux_hpp */
