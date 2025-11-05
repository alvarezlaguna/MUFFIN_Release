//
//  FVM1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 16/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FVM1D_hpp
#define FVM1D_hpp

#include <stdio.h>
#include "SpaceMethod.hpp"
#include "../FluxScheme/FluxScheme.hpp"
#include "../SourceTerm/SourceTerm.hpp"
#include "../SpaceReconstructor/SpaceReconstructor.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"
#include "../Parameters.hpp"
#include <mpi.h>

using namespace Parameters;

class FVM1D : public SpaceMethod
{
public:
    FVM1D(string name) : SpaceMethod(name), m_Fip12(NBEQS), m_Fim12(NBEQS) {}
    ~FVM1D(){}
    void setup();
    void unsetup();
    virtual void setBoundaries();
    virtual void computeRHS();
    virtual double getDtOvDx(){
        int world_size, world_rank;
        // We take the min dtOvdx of all the processors
        double globaldtOvdx;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Allreduce(&m_dtOvdx, &globaldtOvdx, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        return globaldtOvdx;
    };
    virtual double getDt(){
        int world_size, world_rank;
        // We take the min dtOvdx of all the processors
        double globaldt;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Allreduce(&m_dt, &globaldt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        return globaldt;
    };
    
protected:
    unique_ptr<FluxScheme> m_flux;
    unique_ptr<SourceTerm> m_source;
    unique_ptr<SpaceReconstructor> m_reconstructor;
    unique_ptr<BoundaryCondition> m_InletBC;
    unique_ptr<BoundaryCondition> m_OutletBC;
    CellDataRef m_uInlet;
    CellDataRef m_uOutlet;
    vector<double> m_Fip12;
    vector<double> m_Fim12;
};

#endif /* FVM1D_hpp */
