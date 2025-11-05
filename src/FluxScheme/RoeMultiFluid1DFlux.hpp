//
//  RoeMultiFluid1DFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef RoeMultiFluid1DFlux_hpp
#define RoeMultiFluid1DFlux_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "RoeEuler1DFlux.hpp"
#include "../PhysicalModel/EulerEq1D.hpp"
#include "../PhysicalModel/MultiFluid1D.hpp"
#include "HLLEuler1DFlux.hpp"
#include "LFLowMach1DFlux.hpp"


using namespace Parameters;
using namespace std;

class RoeMultiFluid1DFlux : public FluxScheme
{
public:
    RoeMultiFluid1DFlux(string name) : FluxScheme(name), m_flux(NBEQS),
    m_physicalModels(NBFLUIDS),
    m_fluxSchemes(NBFLUIDS),
    m_uRModels(NBFLUIDS),
    m_uLModels(NBFLUIDS),
    m_fluxModels(NBFLUIDS, vector<double>(3, 0.))
    {
        m_fluxSchemes[0].reset(new LFLowMach1DFlux("LFLowMach1DFlux"));
        for(unsigned int iFluid = 1; iFluid < NBFLUIDS; iFluid++) {
            m_fluxSchemes[iFluid].reset(new HLLEuler1DFlux("HLLEuler1DFlux")); // We create a pointer to the fluxes of each
        }
    }
    ~RoeMultiFluid1DFlux(){}
    virtual void setPhysicalModel(PhysicalModel* pm){
        
        FluxScheme::setPhysicalModel(pm);
        MultiFluid1D* model = dynamic_cast<MultiFluid1D*>(m_pm);
        //m_physicalModels = model->getPhysicalModels();
        for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
            PhysicalModel& physicalModel_i = model->getIndividualPhysicalModel(iFluid);
            m_fluxSchemes[iFluid]->setPhysicalModel(&physicalModel_i); // This sets the physical model
        }
    }
    
    vector<double>& operator() (const CellDataRef uL, const CellDataRef uR);
protected:
    vector<double> m_flux;
    vector<unique_ptr<PhysicalModel>> m_physicalModels;                // Vector of pointers to the physical models
    vector<unique_ptr<FluxScheme>> m_fluxSchemes;             // Vector of pointers to the flux Schemes
    vector<CellDataRef> m_uRModels;                       // Vector storing the values of the Right state of the each fluid
    vector<CellDataRef> m_uLModels;                       // Vector storing the values of the Right state of the each fluid
    vector<vector<double>> m_fluxModels;                       // I assume that the models are Euler
};

#endif /* RoeMultiFluid1DFlux_hpp */
