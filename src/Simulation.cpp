//
//  Simulation.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "Simulation.hpp"
#include "MeshData/MeshData.hpp"
#include "MeshData/Cell.hpp"

void Simulation::setup() {

    MeshData& md = MeshData::getInstance();
    md.createData<double>("rhs", NBCELLS*NBEQS);    // Steady residual
    md.createData<double>("physTime", 1);           // Simulation physical time
    
}

void Simulation::unsetup() {
    MeshData& md = MeshData::getInstance();
    //md.createData<Cell>("Cells", NBCELLS);
//    md.deleteData<double>("uCC");
//    md.deleteData<double>("uL");
//    md.deleteData<double>("uR");
    md.deleteData<double>("rhs");
    md.deleteData<double>("physTime");
    
}
