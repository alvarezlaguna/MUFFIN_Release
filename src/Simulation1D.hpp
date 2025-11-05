//
//  Simulate1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Simulation1D_hpp
#define Simulation1D_hpp

#include <stdio.h>
#include "Component.hpp"
#include "Configuration.hpp"
#include "FluxScheme/FluxScheme.hpp"
#include "Simulation.hpp"
#include "SpaceMethod/SpaceMethod.hpp"
#include "TimeMethod/TimeMethod.hpp"
#include "PhysicalModel/PhysicalModel.hpp"
#include "DataWriter/DataWriter.hpp"
#include "CollisionalData/Mixture.hpp"
#include "DiffusiveFlux/DiffusiveFlux.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


namespace py = pybind11;

class Simulation1D : public Simulation {
public:
    Simulation1D(string name) : Simulation(name) {}
    ~Simulation1D(){}
    void setup();                                                       // This function overrides the mother class
    void setupWithOptions(const py::dict &options){
        //import_mpi4py();
        int my_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        
        clock_t start = clock();            // start timing

        if (my_rank == MPI_WRITER){
            py::print( "****************************************************************\n");
            py::print( "\n");
            py::print( "\t\t\tCONFIGURATION PHASE\t\t\n\n");
        }

        Configuration config(options);                 // set the global variables of the simulation
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (my_rank == MPI_WRITER){
            py::print( "****************************************************************\n");
            py::print( "\n");
            py::print( "\t\t\tSET-UP PHASE\t\t\n\n");
        }
        Simulation1D::setup();
    };                               // This function overrides the class up

    void unsetup();                                                     // This function overrides the mother class
    pybind11::array_t<double> simulate();       
    pybind11::array_t<double> advance_one(double dt);                                     // Main FV 1D loop
    void setInitialField();                                             // set initial field
    void writeData(const int iter);                                     // function writing the data to a file
    void computeNorm(vector<double>& norm, vector<double>& normGlobal);                       // function to write the norm
    void addSteps(int nbSteps){NBSTEPS += nbSteps;}                    // function to add steps
    void setResidual(vector<double>& norm){m_Residual = norm;}
    vector<double> getResidual(){return m_Residual;}
    
protected:
    unique_ptr<SpaceMethod> m_spaceMethod;
    unique_ptr<TimeMethod> m_timeMethod;
    unique_ptr<PhysicalModel> m_physicalModel;
    unique_ptr<DataWriter> m_dataWriter;
    unique_ptr<Mixture> m_mixture;
    unique_ptr<DiffusiveFlux> m_diffFlux;
    vector<double> m_Residual;

};

#endif /* Simulation1D_hpp */
