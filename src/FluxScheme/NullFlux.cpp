//
//  NullFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "NullFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
#include <iostream>

using namespace std;

REGISTER_FLUXSCHEME("NullFlux", NullFlux);

vector<double>& NullFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{

    return m_flux;

}


