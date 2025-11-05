//
//  DirichletBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DirichletBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../Parameters.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/MeshData.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <vector>

using namespace Parameters;
using namespace std;

REGISTER_BOUNDARY_CONDITION("Dirichlet", DirichletBC);

CellDataRef DirichletBC::setBoundary(){

    // TODO: Write it general in a class
    for (int iEq = 0; iEq < NBEQS; iEq++){
        if(m_side == "Right"){
            m_value[iEq]   = UOUTLET[iEq];;   // Outlet
        }
        else if(m_side == "Left"){
            m_value[iEq]   = UINLET[iEq];             // Inlet
        }
    }
    return m_value;
}
