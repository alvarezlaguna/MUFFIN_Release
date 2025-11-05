//
//  DiffusiveFluxFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DiffusiveFluxFactory.hpp"

using namespace std;

std::unique_ptr<DiffusiveFlux> DiffusiveFluxFactory::CreateDiffusiveFlux(const string &name)
{
    std::unique_ptr<DiffusiveFlux> DiffusiveFlux;
    if(name == "Null"){DiffusiveFlux.reset(new NullDiffusiveFlux("NullDiffusiveFlux"));}
    if(name == "FourierHeatFlux"){DiffusiveFlux.reset(new FourierHeatFlux("FourierHeatFlux"));}
    else {
        py::print("DiffusiveFlux Model not implemented. Taking No DiffusiveFlux by default\n");
        DiffusiveFlux.reset(new NullDiffusiveFlux("Null"));
    }
    
    return DiffusiveFlux;
}
