//
//  RoeMultiFluidIsothermal1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 02/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef RoeMultiFluidIsothermal1D_hpp
#define RoeMultiFluidIsothermal1D_hpp

#include <stdio.h>
#include <vector>
#include "FluxScheme.hpp"
#include "RoeEulerIsothermal1DFlux.hpp"
#include "HLLEulerIsothermal1DFlux.hpp"
#include "LaxFriedrichFlux.hpp"
#include "../PhysicalModel/EulerIsothermal1D.hpp"
#include "../PhysicalModel/MultiFluidIsothermal1D.hpp"

using namespace Parameters;
using namespace std;

class RoeMultiFluidIsothermal1DFlux : public FluxScheme
{
public:
    RoeMultiFluidIsothermal1DFlux(string name) : FluxScheme(name), m_flux(NBEQS),
    m_physicalModels(NBFLUIDS),
    m_fluxSchemes(NBFLUIDS),
    m_uRModels(NBFLUIDS),
    m_uLModels(NBFLUIDS),
    m_fluxModels(NBFLUIDS, vector<double>(2, 0.))
    {
        for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
            m_fluxSchemes[iFluid].reset(new RoeEulerIsothermal1DFlux("RoeEulerIsothermal1DFlux")); // We create a pointer to the fluxes of each
        }
    }
    ~RoeMultiFluidIsothermal1DFlux(){}
    virtual void setPhysicalModel(PhysicalModel* pm){
        
        FluxScheme::setPhysicalModel(pm);
        MultiFluidIsothermal1D* model = dynamic_cast<MultiFluidIsothermal1D*>(m_pm);
        //m_physicalModels = model->getPhysicalModels();
        for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
            PhysicalModel& physicalModel_i = model->getIndividualPhysicalModel(iFluid);
            m_fluxSchemes[iFluid]->setPhysicalModel(&physicalModel_i); // This sets the physical model and then the speed of sound
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

#endif /* RoeMultiFluidIsothermal1D_hpp */
