//
//  TwoFluidIsothermalSheathBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 03/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermalSheathBC.hpp"
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

REGISTER_BOUNDARY_CONDITION("TwoFluidIsothermalSheath", TwoFluidIsothermalSheathBC);

void TwoFluidIsothermalSheathBC::getWallFlux(const CellDataRef innerU, CellDataRef ghostU){
    //const double massRatio     = MASSRATIO;
    const double soundSpeed_e  = SOUNDSPEED[0]; //We assume the electrons to be the first fluid
    const double pi            = atan(1)*4;
    
    //   Internal values
    const double rhoe_eI        = innerU[0];
    const double u_eI           = innerU[1]/innerU[0];
    //   Ghost values
    const double rhoe_eG        = rhoe_eI;
    //   Wall flux
    const double coeff         = rhoe_eI*soundSpeed_e*sqrt(2/pi);
    const double uOvVth        = abs(u_eI/(sqrt(2)*soundSpeed_e));
    const double rhoeUe_W      = -coeff*(exp(-uOvVth*uOvVth)/2 + sqrt(pi)/2*uOvVth*erfc(-uOvVth));

    //    % Copy values in the boundaries
    ghostU[0] = rhoe_eG;
    ghostU[1] = rhoeUe_W;
    ghostU[2] = innerU[2];
    ghostU[3] = innerU[3];
}

CellDataRef TwoFluidIsothermalSheathBC::setBoundary(){
    vector<Cell1D>& cells      = MeshData::getInstance().getData<Cell1D>("Cells");
    if(m_side == "Right"){// Outlet
        getWallFlux(cells[NBCELLS - 1].uCC, m_value);
        const double ghostFlux_e = -m_value[1]; //We change the sign for the electron flux
        m_value[1]   = 2*ghostFlux_e - cells[NBCELLS - 1].uCC[1]; // We impose the Dirichtlet in the electron flux
    }
    else if(m_side == "Left"){
        getWallFlux(cells[0].uCC, m_value);
        const double ghostFlux_e = m_value[1];
        m_value[1]   = 2*ghostFlux_e - cells[0].uCC[1];  // We impose the Dirichtlet in the electron flux
    }
    
    return m_value;
}
