//
//  TVD2ndOrder3M1DReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVD2ndOrder3M1DReconstructor_hpp
#define TVD2ndOrder3M1DReconstructor_hpp

#include <stdio.h>
#include "TVDReconstructor.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"

class TVD2ndOrder3M1DReconstructor : public TVDReconstructor
{
public:
    TVD2ndOrder3M1DReconstructor(string name): TVDReconstructor(name), m_P_In(NBEQS), m_P_Out(NBEQS), m_P_CC(NBEQS*(NBCELLS)), m_P_L(NBEQS*(NBCELLS)), m_P_R(NBEQS*(NBCELLS)){}
    ~TVD2ndOrder3M1DReconstructor() {}
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
        double p      =  2./3.*((*m_boundaries)[0][2] - 0.5*rho*v*v);

        m_P_In[0] = rho;
        m_P_In[1] = v;
        m_P_In[2] = p;
            
        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            rho    = (*m_cells)[iCell].uCC[0];
            v      = (*m_cells)[iCell].uCC[1]/rho;
            p      = 2./3.*((*m_cells)[iCell].uCC[2] - 0.5*rho*v*v);

            m_P_CC[0*NBCELLS + iCell] = rho;
            m_P_CC[1*NBCELLS + iCell] = v;
            m_P_CC[2*NBCELLS + iCell] = p;
        }
        
        
        rho    =  (*m_boundaries)[1][0];
        v      =  (*m_boundaries)[1][1]/rho;
        p      =  2./3.*((*m_boundaries)[1][2] - 0.5*rho*v*v);

        m_P_Out[0] = rho;
        m_P_Out[1] = v;
        m_P_Out[2] = p;

    }

    void primToCons_LR(){

        for (int iCell = 0; iCell < NBCELLS; iCell++)
        {

            double rho_L    = m_P_L[0*NBCELLS + iCell];
            double v_L      = m_P_L[1*NBCELLS + iCell];
            double p_L      = m_P_L[2*NBCELLS + iCell];


            double rho_R    = m_P_R[0*NBCELLS + iCell];
            double v_R      = m_P_R[1*NBCELLS + iCell];
            double p_R      = m_P_R[2*NBCELLS + iCell];


            (*m_cells)[iCell].uL[0] = rho_L;
            (*m_cells)[iCell].uL[1] = rho_L*v_L;
            (*m_cells)[iCell].uL[2] = 0.5*rho_L*v_L*v_L + 3./2.*p_L;

            (*m_cells)[iCell].uR[0] = rho_R;
            (*m_cells)[iCell].uR[1] = rho_R*v_R;
            (*m_cells)[iCell].uR[2] = 0.5*rho_R*v_R*v_R + 3./2.*p_R;

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

#endif /* TVD2ndOrder3M1DReconstructor_hpp */