//
//  MPIInletBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef MPIInletBC_hpp
#define MPIInletBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class MPIInletBC : public BoundaryCondition
{
public:
    MPIInletBC(string name, string side) : BoundaryCondition(name, side){}
    ~MPIInletBC() {}
    virtual CellDataRef setBoundary();
};

#endif /* MPIInletBC_hpp */
