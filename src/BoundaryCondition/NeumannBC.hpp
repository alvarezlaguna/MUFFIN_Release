//
//  NeumannBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NeumannBC_hpp
#define NeumannBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class NeumannBC : public BoundaryCondition
{
public:
    NeumannBC(string name, string side) : BoundaryCondition(name, side){}
    ~NeumannBC() {}
    virtual CellDataRef setBoundary();
};

#endif /* NeumannBC_hpp */
