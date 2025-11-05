//
//  APEulerFriction1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef APEulerFriction1DFlux_hpp
#define APEulerFriction1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "../PhysicalModel/EulerIsothermal1D.hpp"

using namespace Parameters;
using namespace std;

class APEulerFriction1DFlux : public FluxScheme
{
public:
    APEulerFriction1DFlux(string name) : FluxScheme(name), m_flux(2) {}
    ~APEulerFriction1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};



#endif /* APEulerFriction1DFlux_hpp */
