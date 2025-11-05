//
//  TVD2ndOrderBC1DReconstructor.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 04/12/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TVD2ndOrderBC1DReconstructor.hpp"
#include "../MeshData/Cell1D.hpp"

// Self-register with the factory
#include "ReconstructorRegistrar.hpp"
REGISTER_RECONSTRUCTOR("TVD2ndOrder1DBC", TVD2ndOrderBC1DReconstructor);

void TVD2ndOrderBC1DReconstructor::reconstructField() {
    //get the data to reconstruct
    vector<Cell1D>& cells          = MeshData::getInstance().getData<Cell1D>("Cells");
    // old boundary call  = MeshData::getInstance().getData<double>("boundaries");
    vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");

    //vector<double>& x = MeshData::getInstance().getData<double>("x");
    
    const int n1    = NBCELLS - 1;
    
    double a, b, theta, limiterR, limiterL;
    double uInlet, uOutlet;
    
    // Reference to the limiter (initialized in the TVDReconstructor class)
    Limiter& phi = *m_limiter;
    
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
        
        uInlet  = boundaries[0][iEq];
        uOutlet = boundaries[1][iEq];
        
        a        = cells[0].uCC[iEq] - uInlet;
        b        = cells[1].uCC[iEq] - cells[0].uCC[iEq];
        theta    = (a != 0.0) ? b/a : 0.;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0.;
        limiterL = phi(theta);
        
        //if(uInlet == cells[0].uCC[iEq]) // Case of Neumann boundary
        //{
        //    cells[0].uR[iEq]   = cells[0].uCC[iEq];
        //    cells[0].uL[iEq]   = cells[0].uCC[iEq];
        //}
        //else
        //{
        cells[0].uR[iEq]   = cells[0].uCC[iEq]; // First order in the first and last cell
        cells[0].uL[iEq]   = cells[0].uCC[iEq];
        //}
        
        
        for (int i = 1; i < n1; i++)
        {
            a        = cells[i].uCC[iEq] - cells[i-1].uCC[iEq];
            b        = cells[i+1].uCC[iEq] - cells[i].uCC[iEq];
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            cells[i].uR[iEq]    = cells[i].uCC[iEq] + 0.5*limiterR*(cells[i].uCC[iEq] - cells[i-1].uCC[iEq]); //TODO: Non-homogeneous mesh implementation
            cells[i].uL[iEq]    = cells[i].uCC[iEq] - 0.5*limiterL*(cells[i+1].uCC[iEq] - cells[i].uCC[iEq]);
        }
        a        = cells[n1].uCC[iEq] - cells[n1-1].uCC[iEq];
        b        = uOutlet - cells[n1].uCC[iEq];
        theta    = (a != 0.0) ? b/a : 0;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0;
        limiterL = phi(theta);
        
        //if(uOutlet == cells[n1].uCC[iEq]) // Case of Neumann boundary
        //{
        //    cells[n1].uR[iEq]   = cells[n1].uCC[iEq];
        //    cells[n1].uL[iEq]   = cells[n1].uCC[iEq];
        //}
        //else {
        cells[n1].uR[iEq]   = cells[n1].uCC[iEq]; // First order in the first and last cell
        cells[n1].uL[iEq]   = cells[n1].uCC[iEq];
        //}
    }
}
