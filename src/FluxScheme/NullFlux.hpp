//
//  NullFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NullFlux_hpp
#define NullFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class NullFlux : public FluxScheme
{
public:
    NullFlux(string name) : FluxScheme(name), m_flux(NBEQS) {}
    ~NullFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};

#endif /* NullFlux_hpp */
