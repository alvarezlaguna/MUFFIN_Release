//
//  HLLEulerIsothermalScaled1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef HLLEulerIsothermalScaled1DFlux_hpp
#define HLLEulerIsothermalScaled1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "../PhysicalModel/EulerIsothermal1D.hpp"

using namespace Parameters;
using namespace std;

class HLLEulerIsothermalScaled1DFlux : public FluxScheme
{
public:
    HLLEulerIsothermalScaled1DFlux(string name) : FluxScheme(name), m_flux(2) {}
    ~HLLEulerIsothermalScaled1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};


#endif /* HLLEulerIsothermalScaled1DFlux_hpp */
