
//
//  HLLEuler1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 28/01/20.
//  Copyright © 2020 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef HLLEuler1DFlux_hpp
#define HLLEuler1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class HLLEuler1DFlux : public FluxScheme
{
public:
    HLLEuler1DFlux(string name) : FluxScheme(name), m_flux(3) {}
    ~HLLEuler1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};

#endif /* HLLEuler1DFlux_hpp */
