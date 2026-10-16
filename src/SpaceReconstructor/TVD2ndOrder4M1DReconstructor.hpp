//
//  TVD2ndOrder4M1DReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVD2ndOrder4M1DReconstructor_hpp
#define TVD2ndOrder4M1DReconstructor_hpp

#include <stdio.h>
#include "TVDReconstructor.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"

class TVD2ndOrder4M1DReconstructor : public TVDReconstructor
{
public:
    TVD2ndOrder4M1DReconstructor(string name): TVDReconstructor(name), m_P_In(NBEQS), m_P_Out(NBEQS), m_P_CC(NBEQS*(NBCELLS)), m_P_L(NBEQS*(NBCELLS)), m_P_R(NBEQS*(NBCELLS)){}
    ~TVD2ndOrder4M1DReconstructor() {}
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

        m_P_In[0] = rho;
        m_P_In[1] = v;
        m_P_In[2] = p;
        // m_P_In[2] = p/rho;
        m_P_In[3] = (*m_boundaries)[0][3];
            
        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            rho    = (*m_cells)[iCell].uCC[0];
            v      = (*m_cells)[iCell].uCC[1]/rho;
            p      = (*m_cells)[iCell].uCC[2] - rho*v*v;

            m_P_CC[0*NBCELLS + iCell] = rho;
            m_P_CC[1*NBCELLS + iCell] = v;
            m_P_CC[2*NBCELLS + iCell] = p;
            // m_P_CC[2*NBCELLS + iCell] = p/rho;
            m_P_CC[3*NBCELLS + iCell] = (*m_cells)[iCell].uCC[3];

        }
        
        
        rho    =  (*m_boundaries)[1][0];
        v      =  (*m_boundaries)[1][1]/rho;
        p      =  (*m_boundaries)[1][2] - rho*v*v;

        m_P_Out[0] = rho;
        m_P_Out[1] = v;
        m_P_Out[2] = p;
        // m_P_Out[2] = p/rho;
        m_P_Out[3] = (*m_boundaries)[1][3];

    }

    void primToCons_LR(){

        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            double rho_L    = m_P_L[0*NBCELLS + iCell];
            double v_L      = m_P_L[1*NBCELLS + iCell];
            double p_L      = m_P_L[2*NBCELLS + iCell];
            // double p_L      = m_P_L[2*NBCELLS + iCell]*rho_L;
            double p_perp_L = m_P_L[3*NBCELLS + iCell];

            double rho_R    = m_P_R[0*NBCELLS + iCell];
            double v_R      = m_P_R[1*NBCELLS + iCell];
            double p_R      = m_P_R[2*NBCELLS + iCell];
            // double p_R      = m_P_R[2*NBCELLS + iCell]*rho_R;
            double p_perp_R = m_P_R[3*NBCELLS + iCell];

            (*m_cells)[iCell].uL[0] = rho_L;
            (*m_cells)[iCell].uL[1] = rho_L*v_L;
            (*m_cells)[iCell].uL[2] = rho_L*v_L*v_L + p_L;
            (*m_cells)[iCell].uL[3] = p_perp_L;

            (*m_cells)[iCell].uR[0] = rho_R;
            (*m_cells)[iCell].uR[1] = rho_R*v_R;
            (*m_cells)[iCell].uR[2] = rho_R*v_R*v_R + p_R;
            (*m_cells)[iCell].uR[3] = p_perp_R;

        }


    }

private:
    vector<Cell1D>* m_cells;
    vector<CellDataRef>* m_boundaries ;
    vector<double>* m_x;


    vector<double> m_P_In;
    vector<double> m_P_Out;
    vector<double> m_P_CC;
    vector<double> m_P_L;
    vector<double> m_P_R;
};

#endif /* TVD2ndOrder4M1DReconstructor_hpp */