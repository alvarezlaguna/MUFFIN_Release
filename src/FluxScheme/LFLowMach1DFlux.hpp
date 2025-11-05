//
//  LFLowMach1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 06/03/20.
//  Copyright © 2020 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LFLowMach1DFlux_hpp
#define LFLowMach1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class LFLowMach1DFlux : public FluxScheme
{
public:
    LFLowMach1DFlux(string name) : FluxScheme(name), m_flux(3) {}
    ~LFLowMach1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};
#endif /* LFLowMach1DFlux_hpp */
