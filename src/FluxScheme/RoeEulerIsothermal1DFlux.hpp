//
//  RoeEulerIsothermal1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef RoeEulerIsothermal1DFlux_hpp
#define RoeEulerIsothermal1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "../PhysicalModel/EulerIsothermal1D.hpp"

using namespace Parameters;
using namespace std;

class RoeEulerIsothermal1DFlux : public FluxScheme
{
public:
    RoeEulerIsothermal1DFlux(string name) : FluxScheme(name), m_flux(2) {}
    ~RoeEulerIsothermal1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};

#endif /* RoeEulerIsothermal1DFlux_hpp */
