//
//  main.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/07/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//
#include <iostream>
#include <ctime>
#include <string>

#include "Simulation1D.hpp"
#include "Configuration.hpp"
#include "SourceTerm/InputFilesProvider.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

//#include <mpi4py/mpi4py.h>
#include <petsc.h>

#include <mpi.h>

namespace py = pybind11;
using namespace std;
using namespace Parameters;

// -------------------------------------------------------------------------- //

vector<double> Solver(const py::dict options)
{
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

    Simulation1D sim("Fluid Simulation");  // create simulation object
    sim.setup();                        // generate the mesh data
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (my_rank == MPI_WRITER){
        py::print( "\n****************************************************************\n\n");
        py::print( "\n");
        py::print( "\t\t\tSIMULATION PHASE\t\t\n\n");
        py::print( "\n****************************************************************\n\n");
    }
    
    clock_t start_sim = clock();            // start timing

    sim.simulate();                     // simulate the flow
    clock_t end_sim = clock();   // end timing
    int NBCELLS_Total = 0;
    for (int iP = 0; iP < world_size; ++iP){NBCELLS_Total += NBCELLS_MPI[iP];}
    
    if (my_rank == MPI_WRITER){
        py::print( "Simulation on ",  NBCELLS_Total, " cells took ",
              double(end_sim - start_sim)/CLOCKS_PER_SEC, " s \n");
        py::print( "\n****************************************************************\n\n");
        py::print( "\n");
        py::print( "\t\t\tUNSET-UP PHASE\t\t\n\n");
        py::print( "\n****************************************************************\n\n");
    }

    vector<double> Residual = sim.getResidual();

    sim.unsetup();                      // clean up the mesh data
    
    clock_t end = clock();   // end timing
    if (my_rank == MPI_WRITER){
        py::print( "Simulation on ",  NBCELLS_Total, " cells took ",
                      double(end-start)/CLOCKS_PER_SEC, " s \n");
    }
    // MPI_Finalize();
    return Residual;

}

Simulation1D& Setup(const py::dict options)
{
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

    Simulation1D sim("Fluid Simulation");  // create simulation object
    sim.setup();                        // generate the mesh data
    MPI_Barrier(MPI_COMM_WORLD);

    return sim;
}

PYBIND11_MODULE(muffin, m) 
{
    m.doc() = "pybind11 simulation plugin"; // module docstring
    
    m.def("Solver", &Solver, "CFD simulation of the multi-fluid plasma equations");

    py::class_<Simulation1D>(m, "Simulation1D")
        .def(py::init<string &>())
        .def("simulate", &Simulation1D::simulate,pybind11::return_value_policy::take_ownership)
        .def("advance_one", &Simulation1D::advance_one,pybind11::return_value_policy::take_ownership)
        .def("setup", &Simulation1D::setupWithOptions)
        .def("addSteps", &Simulation1D::addSteps);

    m.def("Setup", &Setup, "Set-up of simulation", pybind11::return_value_policy::reference);

    py::class_<InputFilesProvider>(m, "InputFilesProvider")
        .def_static("getInstance", &InputFilesProvider::getInstance, pybind11::return_value_policy::reference)
        .def("setData", &InputFilesProvider::setData);
}














