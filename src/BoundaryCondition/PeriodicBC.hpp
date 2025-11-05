//
//  PeriodicBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PeriodicBC_hpp
#define PeriodicBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class PeriodicBC : public BoundaryCondition
{
public:
    PeriodicBC(string name, string side) : BoundaryCondition(name, side){}
    ~PeriodicBC() {}
    virtual CellDataRef setBoundary();
};

#endif /* PeriodicBC_hpp */
