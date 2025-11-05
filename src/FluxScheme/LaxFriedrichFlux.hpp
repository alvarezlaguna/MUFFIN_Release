//
//  LaxFriedrichFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LaxFriedrichFlux_hpp
#define LaxFriedrichFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class LaxFriedrichFlux : public FluxScheme
{
public:
    LaxFriedrichFlux(string name) : FluxScheme(name), m_flux(NBEQS), m_flux_l(NBEQS), m_eigenvals_l(NBEQS), m_eigenvals_r(NBEQS), m_flux_r(NBEQS) {}
    ~LaxFriedrichFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
    vector<double> m_eigenvals_l;
    vector<double> m_eigenvals_r;
    vector<double> m_flux_l;
    vector<double> m_flux_r;
};

#endif /* LaxFriedrichFlux_hpp */
