//
//  RoeEuler1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef RoeEuler1DFlux_hpp
#define RoeEuler1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class RoeEuler1DFlux : public FluxScheme
{
public:
    RoeEuler1DFlux(string name) : FluxScheme(name), m_flux(3) {}
    ~RoeEuler1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};


#endif /* RoeEuler1DFlux_hpp */
