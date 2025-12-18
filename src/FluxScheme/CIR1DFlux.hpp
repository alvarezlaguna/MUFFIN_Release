//
//  CIR1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef CIR1DFlux_hpp
#define CIR1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "../PhysicalModel/LinearizedGrad1D.hpp"

using namespace Parameters;
using namespace std;

class CIR1DFlux : public FluxScheme
{
public:
    CIR1DFlux(string name) : FluxScheme(name), m_flux(Parameters::LEVEL) {}
    ~CIR1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
};

#endif /* CIR1DFlux_hpp */
