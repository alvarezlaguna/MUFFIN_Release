//
//  DirichletBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DirichletBC_hpp
#define DirichletBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class DirichletBC : public BoundaryCondition
{
public:
    DirichletBC(string name, string side) : BoundaryCondition(name, side){}
    ~DirichletBC() {}
    virtual CellDataRef setBoundary();
};

#endif /* DirichletBC_hpp */
