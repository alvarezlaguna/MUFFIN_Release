//
//  FVM1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 16/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FVM1D.hpp"
#include "../FluxScheme/FluxSchemeFactory.hpp"
#include "../SourceTerm/SourceTermFactory.hpp"
#include "../SpaceReconstructor/ReconstructorFactory.hpp"
#include "../BoundaryCondition/BoundaryConditionFactory.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"
#include <pybind11/pybind11.h>
#include <mpi.h>


namespace py = pybind11;


void FVM1D::setup(){
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(my_rank == MPI_WRITER)
        py::print("Setting-up Space discretization class:\t ",this->getName(),"\n");
    
    // TODO: Change this

    
    // Initialize the flux scheme
    m_flux = FluxSchemeFactory::CreateFluxScheme(FLUXSCHEMENAME);
    m_flux->setPhysicalModel(m_pm);

    
    // Initialize the source term
    m_source = SourceTermFactory::CreateSourceTerm(SOURCETERMNAME);
    m_source->setPhysicalModel(m_pm);
    m_source->setup();
    
    // Initialize the reconstructor for second order solution in space
    m_reconstructor = ReconstructorFactory::CreateReconstructor(RECONSTRUCTORNAME);
    m_reconstructor->setup(); // initialize the limiter
    
    // TODO: Change this to be set from input file
    m_InletBC  = BoundaryConditionFactory::CreateBoundaryCondition(INLETTYPE,"Left");
    
    m_OutletBC  = BoundaryConditionFactory::CreateBoundaryCondition(OUTLETTYPE,"Right");
    
    if(my_rank == MPI_WRITER){
        py::print("Space discretization using the flux scheme:\t ",m_flux->getName(),"\n");
        py::print("Space discretization using the source term:\t ",m_source->getName(),"\n");
        py::print("Space discretization using the reconstructor:\t ",m_reconstructor->getName(),"\n");
    }
}

void FVM1D::unsetup(){
    m_reconstructor->unsetup();
}

void FVM1D::setBoundaries(){

    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    
    // MPI excahnge of data
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Status status;
    
    int left_neigh = (my_rank>0) ? my_rank-1 : MPI_PROC_NULL;
    int right_neigh = (my_rank<(world_size-1)) ? my_rank+1 : MPI_PROC_NULL;
    

    
    // Send and receive the neighbour information
    vector<double> sendbufL=vector<double>(NBEQS);
    vector<double> sendbufR=vector<double>(NBEQS);
    for(unsigned int iEq = 0; iEq < NBEQS; iEq++){
        sendbufL[iEq] = cells[0].uCC[iEq];
        sendbufR[iEq] = cells[NBCELLS - 1].uCC[iEq];
    }

    MPI_Sendrecv(sendbufL.data(), sendbufL.size(), MPI_DOUBLE, left_neigh, 6, UOUTLET.data(),  UOUTLET.size(), MPI_DOUBLE, right_neigh, 6, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(sendbufR.data(), sendbufR.size(), MPI_DOUBLE, right_neigh, 6, UINLET.data(),  UINLET.size(), MPI_DOUBLE, left_neigh, 6, MPI_COMM_WORLD, &status);
    
    // Here I set-up the m_uInlet
    m_uInlet  = m_InletBC->setBoundary();
    m_uOutlet = m_OutletBC->setBoundary();


// Example of comment to debug the boundaries
//    for (int iP = 0; iP < world_size; ++iP){
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (my_rank == iP) {
//            cout<<"In rank = "<<my_rank<<" boundaries[0] = "<< boundaries[0] <<" boundaries[1] = "<< boundaries[1] <<"\n";
//            cout<<"In rank = "<<my_rank<<" boundaries[2] = "<< boundaries[2] <<" boundaries[3] = "<< boundaries[3] <<"\n";
//            cout<<"In rank = "<<my_rank<<" boundaries[4] = "<< boundaries[4] <<" boundaries[5] = "<< boundaries[5] <<"\n";}
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
}

void FVM1D::computeRHS(){
    
    const int    n1 = NBCELLS - 1;          // number of cells minus 1
    
    FluxScheme& flux1D          = *m_flux;      // reference is used to keep the same syntax as before
    SourceTerm& sourceterm      = *m_source;
    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs         = MeshData::getInstance().getData<double>("rhs");
    
    // copy the boundaries to a value
    //not anymore since it is the same ref
    
    // Reconstruction
    m_reconstructor->reconstructField();
    
    double maxEigenvalue = 0;

    // Main loop
    m_Fip12         = flux1D(cells[0].uR, cells[1].uL);
    m_Fim12         = flux1D(m_uInlet, cells[0].uL);
    // First cell
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS]        = m_Fip12[iEq] - m_Fim12[iEq];

    }
    // compute the max eigenvalue
    double eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[0].uR));
    double eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[1].uL));
    double eigenval_im1R    = abs(m_pm->getMaxEigenvalue(m_uInlet));
    double eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[0].uL));
    maxEigenvalue           = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = CFL/maxEigenvalue;
    m_dt     = m_dtOvdx*cells[0].dx;
    
    // Loop over inner cells
    for (int iCell = 1; iCell < n1; ++iCell) {
        m_Fip12         = flux1D(cells[iCell].uR, cells[iCell+1].uL);
        m_Fim12         = flux1D(cells[iCell-1].uR, cells[iCell].uL);
        for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
            // compute the flux
            rhs[iEq*NBCELLS + iCell]    = m_Fip12[iEq] - m_Fim12[iEq];
        }
        
        // compute the max eigenvalue
        eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[iCell].uR));
        eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[iCell+1].uL));
        eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[iCell-1].uR));
        eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[iCell].uL));
        maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
        m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
        m_dt     = min(m_dt, m_dtOvdx*cells[iCell].dx);
    } // integral flux on internal cells
    // Last cell
    m_Fip12         = flux1D(cells[n1].uR, m_uOutlet);
    m_Fim12         = flux1D(cells[n1-1].uR, cells[n1].uL);
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS + n1]       = m_Fip12[iEq] - m_Fim12[iEq];
    }
    
    // Compute the sourceTerm
    sourceterm.computeSource();
    
    // compute the max eigenvalue
    eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[n1].uR));
    eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(m_uOutlet));
    eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[n1-1].uR));
    eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[n1].uL));
    maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
    m_dt     = min(m_dt, m_dtOvdx*cells[n1].dx);

    const double sourceFrequency = sourceterm.getFrequency();
    if(abs(sourceFrequency) < 1e-12){ // case where the frequency is close to zero
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_convective       = m_dt;
        m_dtOvdx = dtOvdx_convective;
        m_dt     = dt_convective;
    }
    else{
        const double Dx                  = cells[0].dx; //Assuming the cells have the same delta_x
        const double dtOvdx_source       = CFL/(sourceFrequency*Dx);
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_source           = CFL/(sourceFrequency);
        const double dt_convective       = m_dt;

        m_dtOvdx = min(dtOvdx_source, dtOvdx_convective);
        m_dt     = min(dt_source, dt_convective);
    }

}
