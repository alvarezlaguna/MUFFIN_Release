//
//  CIR1DFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "CIR1DFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
using namespace std;
using namespace Parameters;

REGISTER_FLUXSCHEME("CIR1D", CIR1DFlux);

vector<double>& CIR1DFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    const double cols = Parameters::LEVEL;
    for (int i = 0; i < Parameters::LEVEL; ++i) {
        double sum = 0.0;
        const int row_offset = i * cols;  // Cache the offset for the current row
        for (int j = 0; j < cols; ++j) {
            sum += Parameters::A_MATRIX_PLUS_GRAD[row_offset + j] * uL[j] + Parameters::A_MATRIX_MINUS_GRAD[row_offset + j] * uR[j];
        }
        m_flux[i] = sum;
    }
    
    return m_flux;
}