//
//  NoReconstructor.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "NoReconstructor.hpp"
#include "../MeshData/Cell1D.hpp"
// Self-register with the factory
#include "ReconstructorRegistrar.hpp"
REGISTER_RECONSTRUCTOR("1stOrder", NoReconstructor);

void NoReconstructor::reconstructField() {
    //get the data to reconstruct
    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    
    // Copying data of the cell center to the Left and Right states
    for(unsigned int i = 0; i<NBCELLS; i++) {
        cells[i].uL = cells[i].uCC;
        cells[i].uR = cells[i].uCC;
    }
}
