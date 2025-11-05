//
//  MPIOutletBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "MPIOutletBC.hpp"
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

REGISTER_BOUNDARY_CONDITION("MPIOUTLET", MPIOutletBC);

CellDataRef MPIOutletBC::setBoundary(){
    // TODO: Write it general in a class
    for (int iEq = 0; iEq < NBEQS; iEq++){
        m_value[iEq]   = UOUTLET[iEq];
    }
    return m_value;
}