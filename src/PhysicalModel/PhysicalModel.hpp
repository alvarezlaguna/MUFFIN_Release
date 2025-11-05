//
//  PhysicalModel.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PhysicalModel_hpp
#define PhysicalModel_hpp

#include <stdio.h>
#include <vector>
#include "../Component.hpp"
#include "../CollisionalData/Mixture.hpp"
#include "../DiffusiveFlux/DiffusiveFlux.hpp"


using namespace Parameters;
using namespace std;

class PhysicalModel : public Component
{
public :
    PhysicalModel(string name) : Component(name) {}
    virtual ~PhysicalModel(){};
    void setMixture(Mixture* mixture) {m_mixture = mixture;}
    void setDiffusiveFlux(DiffusiveFlux* diffFlux) {m_diffFlux = diffFlux;}
    Mixture* getMixture(){return m_mixture;}
    virtual vector<double>& getEigenvalues(const CellDataRef u) = 0;
    virtual double getMaxEigenvalue(const CellDataRef u) = 0;
    virtual vector<double>& computePhysicalFlux(const CellDataRef u) = 0;
    
protected :
    Mixture* m_mixture;
    DiffusiveFlux* m_diffFlux;
    
    
};

#endif /* PhysicalModel_hpp */
