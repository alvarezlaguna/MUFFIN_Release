//
//  MPIOutletBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef MPIOutletBC_hpp
#define MPIOutletBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class MPIOutletBC : public BoundaryCondition
{
public:
    MPIOutletBC(string name, string side) : BoundaryCondition(name, side){}
    ~MPIOutletBC() {}
    virtual CellDataRef setBoundary();
};


#endif /* MPIOutletBC_hpp */
