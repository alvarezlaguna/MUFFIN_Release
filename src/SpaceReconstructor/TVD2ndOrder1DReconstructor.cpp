//
//  TVD2ndOrder1DReconstructor.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TVD2ndOrder1DReconstructor.hpp"
#include "../MeshData/Cell1D.hpp"

// Self-register with the factory
#include "ReconstructorRegistrar.hpp"
REGISTER_RECONSTRUCTOR("TVD2ndOrder1D", TVD2ndOrder1DReconstructor);

void TVD2ndOrder1DReconstructor::reconstructField() {
    //get the data to reconstruct
    vector<Cell1D>& cells          = MeshData::getInstance().getData<Cell1D>("Cells");
    // old boundary call  = MeshData::getInstance().getData<double>("boundaries");
    vector<CellDataRef>& boundaries  = MeshData::getInstance().getData<CellDataRef>("boundaries");

    vector<double>& x = MeshData::getInstance().getData<double>("x");
    
    const int n1    = NBCELLS - 1;
    
    double a, b, theta, limiterR, limiterL;
    double uInlet, uOutlet, xInlet, xOutlet;
    
    // Reference to the limiter (initialized in the TVDReconstructor class)
    Limiter& phi = *m_limiter;
    
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
    
        uInlet  = boundaries[0][iEq];
        uOutlet = boundaries[1][iEq];
        
        xInlet  = x[0] - cells[0].dx;
        xOutlet = x[n1] + cells[n1].dx;
        
        a        = (cells[0].uCC[iEq] - uInlet)/(x[0] - xInlet);
        b        = (cells[1].uCC[iEq] - cells[0].uCC[iEq])/(x[1] - x[0]);
        theta    = (a != 0.0) ? b/a : 0.;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0.;
        limiterL = phi(theta);
        
        if(uInlet == cells[0].uCC[iEq]) // Case of Neumann boundary
        {
            cells[0].uR[iEq]   = cells[0].uCC[iEq];
            cells[0].uL[iEq]   = cells[0].uCC[iEq];
        }
        else if (INLETTYPE == "Periodic") // We assume that the first and the last cell are of the same size
        {
//            cout<<"Inlet\n";
            a        = cells[0].uCC[iEq] - cells[n1].uCC[iEq];
            b        = cells[1].uCC[iEq] - cells[0].uCC[iEq];
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            cells[0].uR[iEq]    = cells[0].uCC[iEq] + 0.5*limiterR*(cells[0].uCC[iEq] - cells[n1].uCC[iEq]); //TODO: Non-homogeneous mesh implementation
            cells[0].uL[iEq]    = cells[0].uCC[iEq] - 0.5*limiterL*(cells[1].uCC[iEq] - cells[0].uCC[iEq]);
            
        }
        else
        {
            cells[0].uR[iEq]   = cells[0].uCC[iEq] + 0.5*limiterR*(cells[0].uCC[iEq] - uInlet);
            cells[0].uL[iEq]   = cells[0].uCC[iEq] - 0.5*limiterL*(cells[1].uCC[iEq] - cells[0].uCC[iEq]);
        }
        
        
        for (int i = 1; i < n1; i++)
        {
            double x_ip1    = x[i + 1];
            double x_im1    = x[i - 1];
            double x_i      = x[i];
            
            a        = (cells[i].uCC[iEq] - cells[i-1].uCC[iEq])/(x_i - x_im1);
            b        = (cells[i+1].uCC[iEq] - cells[i].uCC[iEq])/(x_ip1 - x_i);
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            cells[i].uR[iEq]    = cells[i].uCC[iEq] + 0.5*limiterR*(cells[i].uCC[iEq] - cells[i-1].uCC[iEq]); //TODO: Non-homogeneous mesh implementation
            cells[i].uL[iEq]    = cells[i].uCC[iEq] - 0.5*limiterL*(cells[i+1].uCC[iEq] - cells[i].uCC[iEq]);
        }
        
        a        = (cells[n1].uCC[iEq] - cells[n1-1].uCC[iEq])/(x[n1] - x[n1 - 1]);
        b        = (uOutlet - cells[n1].uCC[iEq])/(xOutlet - x[n1]);
        theta    = (a != 0.0) ? b/a : 0;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0;
        limiterL = phi(theta);
        
        if(uOutlet == cells[n1].uCC[iEq]) // Case of Neumann boundary
        {
            cells[n1].uR[iEq]   = cells[n1].uCC[iEq];
            cells[n1].uL[iEq]   = cells[n1].uCC[iEq];
        }
        else if (OUTLETTYPE == "Periodic")
        {
//            cout<<"Outlet\n";
            a        = cells[n1].uCC[iEq] - cells[n1-1].uCC[iEq];
            b        = cells[0].uCC[iEq] - cells[n1].uCC[iEq];
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            cells[n1].uR[iEq]    = cells[n1].uCC[iEq] + 0.5*limiterR*(cells[n1].uCC[iEq] - cells[n1 - 1].uCC[iEq]); //TODO: Non-homogeneous mesh implementation
            cells[n1].uL[iEq]    = cells[n1].uCC[iEq] - 0.5*limiterL*(cells[0].uCC[iEq] - cells[n1].uCC[iEq]);
            
        }
        else {
            cells[n1].uR[iEq]   = cells[n1].uCC[iEq] + 0.5*limiterR*(cells[n1].uCC[iEq] - cells[n1-1].uCC[iEq]);
            cells[n1].uL[iEq]   = cells[n1].uCC[iEq] - 0.5*limiterL*(uOutlet - cells[n1].uCC[iEq]);
        }
    }
}
