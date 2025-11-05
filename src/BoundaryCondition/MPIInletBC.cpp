//
//  MPIInletBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "MPIInletBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../Parameters.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/MeshData.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/Cell1D.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <vector>

using namespace Parameters;
using namespace std;

REGISTER_BOUNDARY_CONDITION("MPIINLET", MPIInletBC);

CellDataRef MPIInletBC::setBoundary(){
    
    // TODO: Write it general in a class
    for (int iEq = 0; iEq < NBEQS; iEq++){
        m_value[iEq]   = UINLET[iEq];
    }
    return m_value;
}

