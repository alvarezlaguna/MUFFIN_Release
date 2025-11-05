//
//  HLLEulerIsothermal1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 29/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef HLLEulerIsothermal1DFlux_hpp
#define HLLEulerIsothermal1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "../PhysicalModel/EulerIsothermal1D.hpp"

using namespace Parameters;
using namespace std;

class HLLEulerIsothermal1DFlux : public FluxScheme
{
public:
    HLLEulerIsothermal1DFlux(string name) : FluxScheme(name), m_flux(2) {}
    ~HLLEulerIsothermal1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};

#endif /* HLLEulerIsothermal1DFlux_hpp */
