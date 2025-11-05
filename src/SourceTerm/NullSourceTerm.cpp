//
//  NullSourceTerm.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "NullSourceTerm.hpp"

// Self-register with the factory
#include "SourceTermRegistrar.hpp"
REGISTER_SOURCETERM("NullSourceTerm", NullSourceTerm);
#include "../MeshData/MeshData.hpp"

using namespace std;

void NullSourceTerm::computeSource()
{
    double* source      = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);

    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            source[iEq*NBCELLS+ iCell] = 0.;
        }
    }
}
