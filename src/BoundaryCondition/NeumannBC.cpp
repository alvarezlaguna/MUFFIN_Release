//
//  NeumannBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "NeumannBC.hpp"
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

REGISTER_BOUNDARY_CONDITION("Neumann", NeumannBC);

CellDataRef NeumannBC::setBoundary(){
    vector<Cell1D>& cells      = MeshData::getInstance().getData<Cell1D>("Cells");
    for (int iEq = 0; iEq < NBEQS; iEq++){
        if(m_side == "Right"){
            m_value[iEq]   = cells[NBCELLS - 1].uCC[iEq];   // Outlet
        }
        else if(m_side == "Left"){
            m_value[iEq]   = cells[0].uCC[iEq];             // Inlet
        }
        // DEBUGGING COMMENTS
        //cout<<"m_value["<<iEq<<"] = "<<m_value[iEq]<<"\n";
    }
    return m_value;
}
