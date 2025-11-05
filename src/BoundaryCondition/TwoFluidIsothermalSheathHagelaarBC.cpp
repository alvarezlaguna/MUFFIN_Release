//
//  TwoFluidIsothermalSheathHagelaarBC.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TwoFluidIsothermalSheathHagelaarBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/MeshData.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "../MeshData/Cell1D.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <vector>
#include <math.h>
#include <cmath>

using namespace Parameters;
using namespace std;

REGISTER_BOUNDARY_CONDITION("TwoFluidIsothermalSheathHagelaar", TwoFluidIsothermalSheathHagelaarBC);

void TwoFluidIsothermalSheathHagelaarBC::getWallFlux(const CellDataRef innerU, CellDataRef ghostU){
    //const double massRatio     = MASSRATIO;
    const double soundSpeed_e  = SOUNDSPEED[0]; //We assume the electrons to be the first fluid
//    const double pi            = atan(1)*4;
    
    //   Internal values
    const double rhoe_eI        = innerU[0];
    //   Ghost values
    const double rhoe_eG        = rhoe_eI;
    //   Wall flux
    const double rhoeUe_W      = -rhoe_eG*soundSpeed_e/sqrt(2*PI);
    
    //    % Copy values in the boundaries
    ghostU[0] = rhoe_eG;
    ghostU[1] = rhoeUe_W;
    ghostU[2] = innerU[2];
    ghostU[3] = innerU[3];
}

CellDataRef TwoFluidIsothermalSheathHagelaarBC::setBoundary(){
    vector<Cell1D>& cells      = MeshData::getInstance().getData<Cell1D>("Cells");
    vector<double>& phi        = MeshData::getInstance().getData<double>("Phi");
    
    const double soundSpeed_e  = SOUNDSPEED[0]; //We assume the electrons to be the first fluid
//    const double pi            = atan(1)*4;
    
    double u_Wall = soundSpeed_e/sqrt(2*PI);
    
    if(m_side == "Right"){// Outlet
        getWallFlux(cells[NBCELLS - 1].uCC, m_value);
        
        double u_inner = cells[NBCELLS - 1].uCC[1]/cells[NBCELLS - 1].uCC[0];
        double u_ghost = 2*u_Wall - u_inner;
        
        double density = (u_ghost != 0 && cells[NBCELLS - 1].uCC[1]/u_ghost > 0)? cells[NBCELLS - 1].uCC[1]/u_ghost : 1e-2*cells[NBCELLS - 1].uCC[0];
        
        m_value[0]  = density;
        m_value[1]  = cells[NBCELLS - 1].uCC[1];
        
        // m_value[0]  = cells[NBCELLS - 1].uCC[0];
        // m_value[1]  = 2*-m_value[1] - cells[NBCELLS - 1].uCC[1];
        
        /*
        double n_I = cells[NBCELLS - 1].uCC[0];
        double u_I = cells[NBCELLS - 1].uCC[1]/cells[NBCELLS - 1].uCC[0];
        double phi_I = phi[NBCELLS - 1];

        double u_G   = 2*u_Wall - u_I;
        double phi_wall = PHIOUT;
        double phi_G = 2*phi_wall - phi_I;
//        double n_G   = n_I*exp(-1./MASSRATIO*0.5*(u_G*u_G - u_I*u_I) + (phi_G - phi_I));
        double n_G   = n_I*exp(-1./MASSRATIO*0.5*(u_G*u_G - u_I*u_I) + (phi_G - phi_I) + 1/MASSRATIO*(COLLELECTRONS*SOUNDSPEED[0] + Nu_iz)*u_Wall);

        m_value[0]  = n_G;
        m_value[1]  = n_G*u_G;
        */
    }
    else if(m_side == "Left"){
        getWallFlux(cells[0].uCC, m_value);
        
        double u_inner = cells[0].uCC[1]/cells[0].uCC[0];
        double u_ghost = -2*u_Wall - u_inner;
        
        double density = (u_ghost != 0 && cells[0].uCC[1]/u_ghost > 0)? cells[0].uCC[1]/u_ghost : 1e-2*cells[0].uCC[0];
        
        m_value[0]   = density;
        m_value[1]   = cells[0].uCC[1];
        
        
        // m_value[0]  = cells[0].uCC[0];
        // m_value[1]  = 2*m_value[1] - cells[0].uCC[1];
        
        /*
        double n_I = cells[0].uCC[0];
        double u_I = cells[0].uCC[1]/cells[0].uCC[0];
        double phi_I = phi[0];

        double u_G   = -2*u_Wall - u_I;
	    double phi_wall = PHIIN;
        double phi_G = 2*phi_wall - phi_I;
//        double n_G   = n_I*exp(-1./MASSRATIO*0.5*(u_G*u_G - u_I*u_I) + (phi_G - phi_I));
        double n_G   = n_I*exp(-1./MASSRATIO*0.5*(u_G*u_G - u_I*u_I) + (phi_G - phi_I) - 1/MASSRATIO*(COLLELECTRONS*SOUNDSPEED[0] + Nu_iz)*u_Wall);

        m_value[0]  = n_G;
        m_value[1]  = n_G*u_G;
        */


    }
    
    return m_value;
}
