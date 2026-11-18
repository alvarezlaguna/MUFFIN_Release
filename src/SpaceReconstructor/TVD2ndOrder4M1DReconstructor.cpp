//
//  TVD2ndOrder4M1DReconstructor.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TVD2ndOrder4M1DReconstructor.hpp"
#include "../MeshData/Cell1D.hpp"

// Self-register with the factory
#include "ReconstructorRegistrar.hpp"
REGISTER_RECONSTRUCTOR("TVD2ndOrder4M1D", TVD2ndOrder4M1DReconstructor);


void TVD2ndOrder4M1DReconstructor::reconstructField() {

    
    const int n1    = NBCELLS - 1;
    
    double a, b, theta, limiterR, limiterL;
    double xInlet, xOutlet;
    
    Limiter& phi = *m_limiter;

    TVD2ndOrder4M1DReconstructor::consToPrim_CC();
    
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
            
        xInlet  = (*m_x)[0]  - ((*m_cells).at(0)).dx;
        xOutlet = (*m_x)[n1] + ((*m_cells).at(n1)).dx;
        
        a        = (m_P_CC[iEq*NBCELLS + 0] - m_P_In[iEq])/((*m_x)[0] - xInlet);
        b        = (m_P_CC[iEq*NBCELLS + 1] - m_P_CC[iEq*NBCELLS + 0])/((*m_x)[1] - (*m_x)[0]);
        theta    = (a != 0.0) ? b/a : 0.;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0.;
        limiterL = phi(theta);
        
        if(m_P_In[iEq] == m_P_CC[iEq*NBCELLS + 0])
        {
            m_P_R[iEq*NBCELLS + 0]   = m_P_CC[iEq*NBCELLS + 0];
            m_P_L[iEq*NBCELLS + 0]   = m_P_CC[iEq*NBCELLS + 0];
        }
        else if (INLETTYPE == "Periodic")
        {
            a        = m_P_CC[iEq*NBCELLS + 0] - m_P_CC[iEq*NBCELLS + n1];
            b        = m_P_CC[iEq*NBCELLS + 1] - m_P_CC[iEq*NBCELLS + 0];
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            m_P_R[iEq*NBCELLS + 0]    = m_P_CC[iEq*NBCELLS + 0] + 0.5*limiterR*(m_P_CC[iEq*NBCELLS + 0] - m_P_CC[iEq*NBCELLS + n1]);
            m_P_L[iEq*NBCELLS + 0]    = m_P_CC[iEq*NBCELLS + 0] - 0.5*limiterL*(m_P_CC[iEq*NBCELLS + 1] - m_P_CC[iEq*NBCELLS + 0]);
            
        }
        else
        {
            m_P_R[iEq*NBCELLS + 0]   = m_P_CC[iEq*NBCELLS + 0] + 0.5*limiterR*(m_P_CC[iEq*NBCELLS + 0] - m_P_In[iEq]);
            m_P_L[iEq*NBCELLS + 0]   = m_P_CC[iEq*NBCELLS + 0] - 0.5*limiterL*(m_P_CC[iEq*NBCELLS + 1] - m_P_CC[iEq*NBCELLS + 0]);
        }
        
        
        for (int i = 1; i < n1; i++)
        {
            double x_ip1    = (*m_x)[i + 1];
            double x_im1    = (*m_x)[i - 1];
            double x_i      = (*m_x)[i];
            
            a        = (m_P_CC[iEq*NBCELLS + i] - m_P_CC[iEq*NBCELLS + i - 1])/(x_i - x_im1);
            b        = (m_P_CC[iEq*NBCELLS + i + 1] - m_P_CC[iEq*NBCELLS + i])/(x_ip1 - x_i);
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            m_P_R[iEq*NBCELLS + i ]    = m_P_CC[iEq*NBCELLS + i] + 0.5*limiterR*(m_P_CC[iEq*NBCELLS + i] - m_P_CC[iEq*NBCELLS + i - 1]);
            m_P_L[iEq*NBCELLS + i ]    = m_P_CC[iEq*NBCELLS + i] - 0.5*limiterL*(m_P_CC[iEq*NBCELLS + i + 1] - m_P_CC[iEq*NBCELLS + i]);
        }
        
        a        = (m_P_CC[iEq*NBCELLS + n1] - m_P_CC[iEq*NBCELLS + n1 - 1])/((*m_x)[n1] - (*m_x)[n1 - 1]);
        b        = (m_P_Out[iEq] - m_P_CC[iEq*NBCELLS + n1])/(xOutlet - (*m_x)[n1]);
        theta    = (a != 0.0) ? b/a : 0;
        limiterR = phi(theta);
        theta    = (b != 0.0) ? a/b : 0;
        limiterL = phi(theta);
        
        if(m_P_Out[iEq] == m_P_CC[iEq*NBCELLS + n1])
        {
            m_P_R[iEq*NBCELLS + n1]   = m_P_CC[iEq*NBCELLS + n1];
            m_P_L[iEq*NBCELLS + n1]   = m_P_CC[iEq*NBCELLS + n1];
        }
        else if (OUTLETTYPE == "Periodic")
        {
            a        = m_P_CC[iEq*NBCELLS + n1] - m_P_CC[iEq*NBCELLS + n1 - 1];
            b        = m_P_CC[iEq*NBCELLS + 0] - m_P_CC[iEq*NBCELLS + n1];
            theta    = (a != 0.0) ? b/a : 0;
            limiterR = phi(theta);
            theta    = (b != 0.0) ? a/b : 0;
            limiterL = phi(theta);
            m_P_R[iEq*NBCELLS + n1]    = m_P_CC[iEq*NBCELLS + n1] + 0.5*limiterR*(m_P_CC[iEq*NBCELLS + n1] - m_P_CC[iEq*NBCELLS + n1 - 1]);
            m_P_L[iEq*NBCELLS + n1]    = m_P_CC[iEq*NBCELLS + n1] - 0.5*limiterL*(m_P_CC[iEq*NBCELLS + 0] - m_P_CC[iEq*NBCELLS + n1]);
            
        }
        else {
            m_P_R[iEq*NBCELLS + n1]   = m_P_CC[iEq*NBCELLS + n1] + 0.5*limiterR*(m_P_CC[iEq*NBCELLS + n1] - m_P_CC[iEq*NBCELLS + n1 - 1]);
            m_P_L[iEq*NBCELLS + n1]   = m_P_CC[iEq*NBCELLS + n1] - 0.5*limiterL*(m_P_Out[iEq] - m_P_CC[iEq*NBCELLS + n1]);
        }
    }

    TVD2ndOrder4M1DReconstructor::primToCons_LR();
}