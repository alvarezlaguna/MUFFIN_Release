//
//  SixMomentAQMOM1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SixMomentAQMOM1DFlux_hpp
#define SixMomentAQMOM1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"

using namespace Parameters;
using namespace std;

class SixMomentAQMOM1DFlux : public FluxScheme
{
public:
    SixMomentAQMOM1DFlux(string name) : FluxScheme(name), m_flux(6), m_flux_l(6), m_eigenvals_l(6), m_eigenvals_r(6), m_flux_r(6) {}
    ~SixMomentAQMOM1DFlux(){}
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
    vector<double> m_eigenvals_l;
    vector<double> m_eigenvals_r;
    vector<double> m_flux_l;
    vector<double> m_flux_r;
};

#endif /* SixMomentAQMOM1DFlux_hpp */
