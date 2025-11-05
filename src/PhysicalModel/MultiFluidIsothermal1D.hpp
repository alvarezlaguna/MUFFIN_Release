//
//  MultiFluidIsothermal1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef MultiFluidIsothermal1D_hpp
#define MultiFluidIsothermal1D_hpp

#include <stdio.h>
#include <iostream>
#include "PhysicalModel.hpp"
#include "../Parameters.hpp"
#include "EulerIsothermal1D.hpp"

using namespace std;
using namespace Parameters;

class MultiFluidIsothermal1D : public PhysicalModel
{
public:
    MultiFluidIsothermal1D(string name) : PhysicalModel(name), m_physFlux(NBEQS),
    m_eigenvals(NBEQS), m_physicalModels(NBFLUIDS), m_uModels(NBFLUIDS),
    m_eigenModels(NBFLUIDS, vector<double>(2, 0.)),
    m_fluxModels(NBFLUIDS, vector<double>(2, 0.))
    {
        for(unsigned int iFluid = 0; iFluid < NBFLUIDS; iFluid++) {
            m_physicalModels[iFluid].reset(new EulerIsothermal1D("EulerIsothermal1D", SOUNDSPEED[iFluid])); // We create a pointer to the physical models of the fluids with a different sound speed.
        }
    }
    ~MultiFluidIsothermal1D(){m_physFlux.clear();}
    vector<unique_ptr<PhysicalModel>> const& getPhysicalModels() const {return m_physicalModels;}
    virtual vector<double>& getEigenvalues(const CellDataRef u);
    virtual double getMaxEigenvalue(const CellDataRef u);
    virtual vector<double>& computePhysicalFlux(const CellDataRef u);
    PhysicalModel& getIndividualPhysicalModel(int iFluid){ return *m_physicalModels[iFluid]; }
    
protected:
    vector<double> m_physFlux;
    vector<double> m_eigenvals;
    vector<unique_ptr<PhysicalModel>> m_physicalModels;       // Vector of pointers to the physical models
    vector<CellDataRef> m_uModels;                       // I assume that the models are Euler
    vector<vector<double>> m_eigenModels;                   // Copy of the eigenvalues to not compute them again. I assume that the models are Euler
    vector<vector<double>> m_fluxModels;                   // Copy of the eigenvalues to not compute them again. I assume that the models are Euler
    
    
};

#endif /* MultiFluidIsothermal1D_hpp */
