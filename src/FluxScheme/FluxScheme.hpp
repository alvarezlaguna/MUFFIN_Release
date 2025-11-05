//
//  FluxScheme.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FluxScheme_hpp
#define FluxScheme_hpp

#include <stdio.h>
#include "../Component.hpp"
#include "../PhysicalModel/PhysicalModel.hpp"

class FluxScheme : public Component
{
public:
    FluxScheme(string name) : Component(name) {}
    virtual ~FluxScheme() {}
    virtual void setPhysicalModel(PhysicalModel* pm){m_pm = pm;}
    virtual vector<double>& operator()(const CellDataRef uL, const CellDataRef uR) = 0;
protected:
    PhysicalModel* m_pm;
};

#endif /* FluxScheme_hpp */
