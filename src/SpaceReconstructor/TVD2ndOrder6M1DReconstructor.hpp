//
//  TVD2ndOrder6M1DReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVD2ndOrder6M1DReconstructor_hpp
#define TVD2ndOrder6M1DReconstructor_hpp

#include <stdio.h>
#include "TVDReconstructor.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"

class TVD2ndOrder6M1DReconstructor : public TVDReconstructor
{
public:
    TVD2ndOrder6M1DReconstructor(string name): TVDReconstructor(name), m_P_In(NBEQS), m_P_Out(NBEQS), m_P_CC(NBEQS*(NBCELLS)), m_P_L(NBEQS*(NBCELLS)), m_P_R(NBEQS*(NBCELLS)){}
    ~TVD2ndOrder6M1DReconstructor() {}
    void setup(){
        TVDReconstructor::setup();
        m_cells    = &MeshData::getInstance().getData<Cell1D>("Cells");
        m_boundaries = &MeshData::getInstance().getData<CellDataRef>("boundaries");
        m_x = &MeshData::getInstance().getData<double>("x");
    };
    void unsetup(){TVDReconstructor::unsetup();}

    virtual void reconstructField();
    void consToPrim_CC(){
            
        double rho    =  (*m_boundaries)[0][0];
        double v      =  (*m_boundaries)[0][1]/rho;
        double p      =  (*m_boundaries)[0][2] - rho*v*v;
        double q      =  (*m_boundaries)[0][3] - rho*v*v*v - 3*v*p;
        double r      =  (*m_boundaries)[0][4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

        m_P_In[0] = rho;
        m_P_In[1] = v;
        m_P_In[2] = p;
        m_P_In[3] = q;
        m_P_In[4] = r;        
        m_P_In[5] = (*m_boundaries)[0][5];        
            
        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            rho    = (*m_cells)[iCell].uCC[0];
            v      = (*m_cells)[iCell].uCC[1]/rho;
            p      = (*m_cells)[iCell].uCC[2] - rho*v*v;
            q      = (*m_cells)[iCell].uCC[3] - rho*v*v*v - 3*v*p;
            r      = (*m_cells)[iCell].uCC[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

            m_P_CC[0*NBCELLS + iCell] = rho;
            m_P_CC[1*NBCELLS + iCell] = v;
            m_P_CC[2*NBCELLS + iCell] = p;
            m_P_CC[3*NBCELLS + iCell] = q;
            m_P_CC[4*NBCELLS + iCell] = r;
            m_P_CC[5*NBCELLS + iCell] = (*m_cells)[iCell].uCC[5];

        }
        
        
        rho    =  (*m_boundaries)[1][0];
        v      =  (*m_boundaries)[1][1]/rho;
        p      =  (*m_boundaries)[1][2] - rho*v*v;
        q      =  (*m_boundaries)[1][3] - rho*v*v*v - 3*v*p;
        r      =  (*m_boundaries)[1][4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

        m_P_Out[0] = rho;
        m_P_Out[1] = v;
        m_P_Out[2] = p;
        m_P_Out[3] = q;
        m_P_Out[4] = r;   
        m_P_Out[5] = (*m_boundaries)[1][5];

    }

    void primToCons_LR(){

        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            double rho_L    = m_P_L[0*NBCELLS + iCell];
            double v_L      = m_P_L[1*NBCELLS + iCell];
            double p_L      = m_P_L[2*NBCELLS + iCell];
            double q_L      = m_P_L[3*NBCELLS + iCell];
            double r_L      = m_P_L[4*NBCELLS + iCell];
            double p_perp_L = m_P_L[5*NBCELLS + iCell];

            double rho_R    = m_P_R[0*NBCELLS + iCell];
            double v_R      = m_P_R[1*NBCELLS + iCell];
            double p_R      = m_P_R[2*NBCELLS + iCell];
            double q_R      = m_P_R[3*NBCELLS + iCell];
            double r_R      = m_P_R[4*NBCELLS + iCell];
            double p_perp_R = m_P_R[5*NBCELLS + iCell];

            (*m_cells)[iCell].uL[0] = rho_L;
            (*m_cells)[iCell].uL[1] = rho_L*v_L;
            (*m_cells)[iCell].uL[2] = rho_L*v_L*v_L + p_L;
            (*m_cells)[iCell].uL[3] = rho_L*v_L*v_L*v_L + 3*v_L*p_L + q_L;
            (*m_cells)[iCell].uL[4] = rho_L*v_L*v_L*v_L*v_L + 6*v_L*v_L*p_L + 4*v_L*q_L + r_L;
            (*m_cells)[iCell].uL[5] = p_perp_L;

            (*m_cells)[iCell].uR[0] = rho_R;
            (*m_cells)[iCell].uR[1] = rho_R*v_R;
            (*m_cells)[iCell].uR[2] = rho_R*v_R*v_R + p_R;
            (*m_cells)[iCell].uR[3] = rho_R*v_R*v_R*v_R + 3*v_R*p_R + q_R;
            (*m_cells)[iCell].uR[4] = rho_R*v_R*v_R*v_R*v_R + 6*v_R*v_R*p_R + 4*v_R*q_R + r_R;
            (*m_cells)[iCell].uR[5] = p_perp_R;

        }


    }

// void primToCons(){
            
//         double rho    =  m_boundaries[0][0];
//         double v      =  m_boundaries[0][1]/rho;
//         double p      =  m_boundaries[0][2] - rho*v*v;
//         double q      =  m_boundaries[0][3] - rho*v*v*v - 3*v*p;
//         double r      =  m_boundaries[0][4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

//         m_P_In[0] = rho;
//         m_P_In[1] = v;
//         m_P_In[2] = p;
//         m_P_In[3] = q;
//         m_P_In[4] = r;        
            
//         for (int iCell = 0; i < NBCELLS; iCell++)
//         {

//             rho    = cells[iCell].uCC[0];
//             v      = cells[iCell].uCC[1]/rho;
//             p      = cells[iCell].uCC[2] - rho*v*v;
//             q      = cells[iCell].uCC[3] - rho*v*v*v - 3*v*p;
//             r      = cells[iCell].uCC[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

//             m_P_CC[0*NBCELLS + iCell] = rho;
//             m_P_CC[1*NBCELLS + iCell] = v;
//             m_P_CC[2*NBCELLS + iCell] = p;
//             m_P_CC[3*NBCELLS + iCell] = q;
//             m_P_CC[4*NBCELLS + iCell] = r;

//         }
        
        
//         rho    =  m_boundaries[1][0];
//         v      =  m_boundaries[1][1]/rho;
//         p      =  m_boundaries[1][2] - rho*v*v;
//         q      =  m_boundaries[1][3] - rho*v*v*v - 3*v*p;
//         r      =  m_boundaries[1][4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;

//         m_P_Out[0] = rho;
//         m_P_Out[1] = v;
//         m_P_Out[2] = p;
//         m_P_Out[3] = q;
//         m_P_Out[4] = r;     

    // }

private:
    vector<Cell1D>* m_cells; // Pointer to the cells to avoid lookup
    vector<CellDataRef>* m_boundaries ;
    vector<double>* m_x;


    vector<double> m_P_In;
    vector<double> m_P_Out;
    vector<double> m_P_CC;
    vector<double> m_P_L;
    vector<double> m_P_R;
};

#endif /* TVD2ndOrder6M1DReconstructor_hpp */
