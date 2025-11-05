//
//  TwoFluidIsothermalSheathAPHagelaarBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 29/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermalSheathAPHagelaarBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/MeshData.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/Cell1D.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <vector>
#include <math.h>
#include <cmath>

using namespace Parameters;
using namespace std;

REGISTER_BOUNDARY_CONDITION("TwoFluidIsothermalSheathAPHagelaar", TwoFluidIsothermalSheathAPHagelaarBC);

void TwoFluidIsothermalSheathAPHagelaarBC::getWallFlux(const CellDataRef innerU, CellDataRef ghostU){
    //const double massRatio     = MASSRATIO;
    const double soundSpeed_e  = SOUNDSPEED[0]; //We assume the electrons to be the first fluid
    const double pi            = atan(1)*4;
    
    //   Internal values
    const double rhoe_eI        = innerU[0];
    //   Ghost values
    const double rhoe_eG        = rhoe_eI;
    //   Wall flux
    const double rhoeUe_W      = -rhoe_eG*soundSpeed_e/sqrt(2*pi);
    
    //    % Copy values in the boundaries
    ghostU[0] = rhoe_eG;
    ghostU[1] = rhoeUe_W;
    ghostU[2] = innerU[2];
    ghostU[3] = innerU[3];
}

CellDataRef TwoFluidIsothermalSheathAPHagelaarBC::setBoundary(){
    vector<Cell1D>& cells      = MeshData::getInstance().getData<Cell1D>("Cells");
    if(m_side == "Right"){// Outlet
        getWallFlux(cells[NBCELLS - 1].uCC, m_value);
        const double ghostFlux_e = -m_value[1]; //We change the sign for the electron flux
        //m_value[0]   = 2*cells[NBCELLS - 1].uCC[2]*1. - cells[NBCELLS - 1].uCC[0];
        m_value[1]   = 2*cells[NBCELLS - 1].uCC[0]*1. - cells[NBCELLS - 1].uCC[1];
        m_value[2]   = 2*cells[NBCELLS - 1].uCC[0]*1. - cells[NBCELLS - 1].uCC[1];
        //m_value[3]   = 2*cells[NBCELLS - 1].uCC[2]*1. - cells[NBCELLS - 1].uCC[3]; // We impose the Dirichtlet in the electron flux
    }
    else if(m_side == "Left"){
        getWallFlux(cells[0].uCC, m_value);
        const double ghostFlux_e = m_value[1];
        m_value[0]   = 2*cells[0].uCC[0]*1. - cells[0].uCC[0];
        m_value[1]   = 2*cells[0].uCC[1]*1. - cells[0].uCC[1];
        //m_value[3]   = -2*cells[0].uCC[2]*1. - cells[0].uCC[3]; // We impose the Dirichtlet in the electron flux
    }
    
    return m_value;
}
