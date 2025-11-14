//
//  Configuration.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "Configuration.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "CollisionalData/CollisionalData.hpp"
#include "Options.hpp"

//#include <mpi4py/mpi4py.h>
#include <petsc.h>

#include <mpi.h>

namespace py = pybind11;
using namespace std;
using namespace Parameters;

// -------------------------------------------------------------------------- //
Configuration::Configuration(const py::dict &options)
{
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    //import_mpi4py();
    //py::object py_comm = options["MPI_COMM"];
    //PyObject* py_obj = py_comm.ptr();
    //int size, rank;
    //MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
    //MPI_Comm_size(&comm_p, &size);
    //MPI_Comm_rank(&comm_p, &rank);
    
    //py::print("Hello from ", world_rank);
    
    //PetscInitializeNoArguments();
    //PetscMPIInt    rank, size;
    //MPI_Comm_size(PETSC_COMM_WORLD,&size);
    //MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    //PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n", size, rank);
    
    
    // Converting types
    // TODO: check the compatibility of the conditions
    
// Loop to check all the options
//    for (auto item : options)
//    {
//        std::cout << "key: " << item.first << ", value=" << item.second << std::endl;
//
//
//    };
    // Read the mixture
    py::object getField  = options.attr("get");
    // Check if we have the mixture option
    if(!options.attr("get")("mixture").is_none()){
        auto mixture     = options.attr("get")("mixture"); //handle with the python mixture
        auto mixData_py  = mixture.attr("get")("mixtureData");    //handle with the python mixture data
        auto elecData_py = mixture.attr("get")("electronData");   //handle with the data
        auto ionData_py  = mixture.attr("get")("ionData");       //handle with the data
        auto p_gas_py   = mixture.attr("get")("p_gas");       //handle with the data
        auto T_gas_py  = mixture.attr("get")("T_gas");       //handle with the data
        // Copy the mixture in the collisional Data singleton
        CollisionalData& cd = CollisionalData::getInstance();
        cd.createData<string>("mixtureType");        //Create the name in memory
        cd.createData<py::handle>("mixtureData");           //Create the name in memory
        cd.createData<py::handle>("electronData");   //Create the data in memory
        cd.createData<py::handle>("ionData");        //Create the data in memory
        cd.createData<py::handle>("p_gas");        //Create the data in memory
        cd.createData<py::handle>("T_gas");        //Create the data in memory
        auto& mixType  = CollisionalData::getInstance().getData<string>("mixtureType");
        auto& mixData  = CollisionalData::getInstance().getData<py::handle>("mixtureData");
        auto& elecData = CollisionalData::getInstance().getData<py::handle>("electronData");
        auto& ionData  = CollisionalData::getInstance().getData<py::handle>("ionData");
        auto& p_gas   = CollisionalData::getInstance().getData<py::handle>("p_gas");
        auto& T_gas   = CollisionalData::getInstance().getData<py::handle>("T_gas");
        mixData  = mixData_py;  // copy the python mixture object
        elecData = elecData_py; // copy the python list of objects
        ionData  = ionData_py;  // copy the python list of objects
        p_gas    = p_gas_py;
        T_gas    = T_gas_py;
        mixType  = py::str(mixture.attr("get")("type"));
    }
    if(options.attr("get")("mixture").is_none()){ // No mixture defined
        CollisionalData& cd = CollisionalData::getInstance();
        cd.createData<string>("mixtureType");        //Create the name in memory
        auto& mixType  = CollisionalData::getInstance().getData<string>("mixtureType");
        mixType        = "Null";
    }
    
    // Check if we have the diffusive flux option
    if(!options.attr("get")("diffusiveFlux").is_none()){
        auto diffFlux    = options.attr("get")("diffusiveFlux");

        // Copy the diffusion model in the options Data singleton
        Options& op = Options::getInstance();
        op.createData<string>("DiffusionType");        //Create the name in memory
        auto& diffModel = Options::getInstance().getData<string>("DiffusionType");

        diffModel  = py::str(diffFlux.attr("get")("type"));
    }
    else if(options.attr("get")("diffusiveFlux").is_none()){ // No diffusion model defined
        Options& op = Options::getInstance();
        op.createData<string>("DiffusionType");        //Create the name in memory
        auto& diffModel = Options::getInstance().getData<string>("DiffusionType");
        diffModel        = "Null";
    }
    
    //py::object getField = options.attr("get");
    py::int_ nbEqs_py = getField("nbEqs");
    int nbEqs = nbEqs_py;
    
    py::int_ nbFluids_py = getField("nbFluids");
    int nbFluids = nbFluids_py;

    
    py::list nbCells_py = options["nbCells"];
    py::int_ nbCells_i_py = nbCells_py[world_rank];
    int nbCells = nbCells_i_py;

    
    Parameters::NBCELLS_MPI.resize(world_size);
    for (int iP = 0; iP < world_size; ++iP){
        py::int_ nbCells_i_py = nbCells_py[iP];
        Parameters::NBCELLS_MPI[iP] = nbCells_i_py;
    }
    
    py::str geometry_py = options["geometry"];
    string geometry = geometry_py;
    if(geometry == "1D"){
        Parameters::MESH.resize(nbCells + 1);
        py::list mesh_py   = options["mesh"];
        py::float_ length_py = options["length"];
        Parameters::LENGTH = length_py;
        for (int i = 0; i < nbCells + 1; i++){
            py::float_ mesh_py_i = mesh_py[i];
            Parameters::MESH[i] = mesh_py_i;
        }
    }

    
    py::float_ CFL_py = options["CFL"];
    double Courant = CFL_py;

    
    py::dict stopCondition = options["stopCondition"];
    py::str stopCondType_py = stopCondition["type"];
    string stopCondType = stopCondType_py;
    if(stopCondType == "Residual") {
        py::float_ EPS_py = stopCondition["value"];
        Parameters::EPS   = EPS_py;                  // residual threshold for convergence
        
    }

    if(stopCondType == "nbSteps") {
        py::int_ NBSTEPS_py     = stopCondition["value"];
        Parameters::NBSTEPS     = NBSTEPS_py;                  // Number of steps
    }

    if(stopCondType == "Time") {
        py::float_ finalTime_Py     = stopCondition["value"];
        Parameters::FINAL_TIME     = finalTime_Py;                  // Number of steps
    }
    
    // Parameters::DATAWRITER      = "DataWriter1DMPI";
    Parameters::DATAWRITER      = "DataWriter1DH5Py";
    
    py::int_ saveRate_py = options["saveRate"];
    int saveRate = saveRate_py;

    
    py::str resultDir_py = options["resultDir"];
    string resultDir = resultDir_py;
    
    py::str LimiterName_py = options["limiter"];
    string LimiterName = LimiterName_py;
    
    py::str ReconstructorName_py = options["reconstruction"];
    string ReconstructorName = ReconstructorName_py;
    
    py::str FluxSchemeName_py = options["fluxScheme"];
    string FluxSchemeName = FluxSchemeName_py;
    
    py::str TimeSchemeName_py = options["timeScheme"];
    string TimeSchemeName = TimeSchemeName_py;
    
    // Physical Model options
    py::dict PhysicalModel = options["PhysicalModel"];
    py::str  PhysicalModelName_py = PhysicalModel["type"];
    string PhysicalModelName = PhysicalModelName_py;
    if (PhysicalModelName == "AdvectionEq1D") {
        py::float_ A_py = PhysicalModel["A"];
        double speed = A_py;
        Parameters::A = speed; // advection speed of the model
    }
    else if (PhysicalModelName == "EulerEq1D" ) {
        py::float_ gamma_py = PhysicalModel["gamma"];
        double gamma = gamma_py;
        Parameters::GAMMA = gamma; // advection speed of the model
    }
    else if (PhysicalModelName == "MultiFluid1D"){
        Parameters::POLYTROPICINDEX.resize(nbFluids);
        py::list polytropicIndex_py = PhysicalModel["gamma"];
        for (int iFluid = 0; iFluid < nbFluids; ++iFluid){
            py::float_ polytropicIndex_py_i = polytropicIndex_py[iFluid];
            Parameters::POLYTROPICINDEX[iFluid] = polytropicIndex_py_i;
        }
    }

    else if (PhysicalModelName == "EulerIsothermal1D" 
        || PhysicalModelName == "MultiFluidIsothermal1D") {
        Parameters::SOUNDSPEED.resize(nbFluids);
        py::list soundSpeed_py = PhysicalModel["soundSpeed"];
        for (int iFluid = 0; iFluid < nbFluids; ++iFluid){
            py::float_ soundSpeed_py_i = soundSpeed_py[iFluid];
            Parameters::SOUNDSPEED[iFluid] = soundSpeed_py_i;
        }
    }

    // Additional PhysicalModel-specific callbacks (Python-based models)
    if(PhysicalModelName == "PythonPhysicalModel"){
        CollisionalData& cd = CollisionalData::getInstance();
        // Flux function
        py::function flux_function_py = PhysicalModel["flux_function"];
        cd.createData<py::function>("PM_fluxFunction");        //Create the name in memory
        auto& flux_function  = CollisionalData::getInstance().getData<py::function>("PM_fluxFunction");
        flux_function = flux_function_py; 
        // Max eigenvalue function
        py::function MaxEigen_function_py = PhysicalModel["maxEigen_function"];
        cd.createData<py::function>("PM_maxEigen_function");        //Create the name in memory
        auto& maxEigen_function  = CollisionalData::getInstance().getData<py::function>("PM_maxEigen_function");
        maxEigen_function = MaxEigen_function_py; 
    }

    
    py::dict SourceTerm = options["SourceTerm"];
    py::str SourceTermName_py = SourceTerm["type"];
    string SourceTermName = SourceTermName_py;
    if(SourceTermName == "TwoFluidIsothermal1D" 
        || SourceTermName == "TwoFluidIsothermalMPI1D" 
        || SourceTermName == "TwoFluidIsothermal1DPeriodic" )
    {
        py::float_ massRatio_py     = SourceTerm["massRatio"];
        double massRatio            = massRatio_py;
        py::float_ DebyeLength_py   = SourceTerm["DebyeLength"];
        double DebyeLength          = DebyeLength_py;
        py::float_ ionizConst_py    = SourceTerm["ionizConst"];
        double ionizConst = ionizConst_py;
        py::float_ epsIoniz_py      = SourceTerm["epsIoniz"];
        double epsIoniz = epsIoniz_py;
        py::float_ CollIons_py      = SourceTerm["CollIons"];
        double CollIons             = CollIons_py;
        py::float_ CollElectrons_py = SourceTerm["CollElectrons"];
        double CollElectrons = CollElectrons_py;
        py::float_ PhiIn_py         = SourceTerm["PhiIn"];
        double PhiIn                = PhiIn_py;
        py::float_ PhiOut_py        = SourceTerm["PhiOut"];
        double PhiOut               = PhiOut_py;
        
        Parameters::PHIINITIAL.resize(nbCells);
        py::list PhyInitial_py = SourceTerm["PhiInitial"];
        for (int i = 0; i < nbCells; i++){
            py::float_ PhyInitial_py_i = PhyInitial_py[i];
            Parameters::PHIINITIAL[i] = PhyInitial_py_i;
        }
        
        Parameters::MASSRATIO       = massRatio;
        Parameters::DEBYELENGTH     = DebyeLength;
        Parameters::IONIZCONST      = ionizConst;
        Parameters::EPSIONIZ        = epsIoniz;
        Parameters::COLLIONS        = CollIons;
        Parameters::COLLELECTRONS   = CollElectrons;
        Parameters::PHIIN           = PhiIn;
        Parameters::PHIOUT          = PhiOut;
        Parameters::DATAWRITER      = "DataWriter1DH5Py";
    }
    
    Parameters::UINLET.resize(nbEqs);
    // Only rank = 0 has the physical inlet
    if (world_rank == 0){
        py::dict Inlet = options["Inlet"];
        py::str inletType_py = Inlet["type"];
        string inletType = inletType_py;
        if(inletType == "Dirichlet"){
            py::list uInlet_py = Inlet["value"];
            for (int iEq = 0;iEq < nbEqs;++iEq){
                py::float_ uInlet_i = uInlet_py[iEq];
                Parameters::UINLET[iEq] = uInlet_i;           // solution in the inlet
            }
        }
        if (inletType == "PythonBC"){
            CollisionalData& cd = CollisionalData::getInstance();
            auto function_py = Inlet["Function"];
            cd.createData<py::function>("functionInlet");        //Create the name in memory
            auto& function  = CollisionalData::getInstance().getData<py::function>("functionInlet");
            
            function = function_py;

        }
        if (inletType == "WallIodineTemperatureBC"){
            py::float_ gamma_wall_BC_py  = Inlet["gamma_wall_BC"];
            double gamma_wall_BC         = gamma_wall_BC_py;
            Parameters::GAMMA_WALL_BC   = gamma_wall_BC;
        }
        Parameters::INLETTYPE= inletType;                     // solution in the inlet
    }
    else{
        Parameters::INLETTYPE = "MPIINLET";                   // solution in the inlet
    }
    
    Parameters::UOUTLET.resize(nbEqs);
    // Only rank = last has the physical inlet
    if (world_rank == world_size - 1){
        py::dict Outlet = options["Outlet"];
        py::str outletType_py = Outlet["type"];
        string outletType = outletType_py;
        if(outletType == "Dirichlet"){
            py::list uOutlet_py = Outlet["value"];
            for (int iEq = 0;iEq < nbEqs;++iEq){
                py::float_ uOutlet_i = uOutlet_py[iEq];
                Parameters::UOUTLET[iEq] = uOutlet_i;         // solution in the outlet
            }
        }        
        if (outletType == "PythonBC"){
            CollisionalData& cd = CollisionalData::getInstance();
            auto function_py = Outlet["Function"];
            cd.createData<py::function>("functionOutlet");        //Create the name in memory
            auto& function  = CollisionalData::getInstance().getData<py::function>("functionOutlet");
            function = function_py;

        }
        
        Parameters::OUTLETTYPE= outletType;                   // solution in the outlet
    }
    else{
        Parameters::OUTLETTYPE = "MPIOUTLET";                   // solution in the outlet
    }
    
   
    if(SourceTermName == "PythonSourceTerm"){
        CollisionalData& cd = CollisionalData::getInstance();
        py::function function_py = SourceTerm["Function"];
        cd.createData<py::function>("functionST");        //Create the name in memory
        auto& function  = CollisionalData::getInstance().getData<py::function>("functionST");
        function = function_py; 
    }

    // Setting additional parameters (TODO: refactor)

    if(!options.attr("get")("TGAS").is_none()){
            py::float_ tgas     = options.attr("get")("TGAS");
            Parameters::TGAS = tgas;
    }
    if(!options.attr("get")("TWALL").is_none()){
            py::float_ twall     = options.attr("get")("TWALL");
            Parameters::TWALL = twall;
    }
    if(!options.attr("get")("PGAS").is_none()){
            py::float_ pgas     = options.attr("get")("PGAS");
            Parameters::PGAS = pgas;
    }

    // General options

    py::list InitialField_py = options["initialField"];
    
    // Changing the values in the namespace that is used in the simulation
    Parameters::NBEQS               = nbEqs;                    // number of equations
    Parameters::NBCELLS             = nbCells;                  // number of cells
    Parameters::NBFLUIDS            = nbFluids;                 // number of fluids
    Parameters::CFL                 = Courant;                  // Courant Friedrichs Lewy number
    Parameters::SAVERATE            = saveRate;                 // save rate
    Parameters::RESULTDIRECTORY     = resultDir;                // result directory
    Parameters::LIMITERNAME         = LimiterName;              // Limiter Name
    Parameters::RECONSTRUCTORNAME   = ReconstructorName;        // Reconstructor Name
    Parameters::FLUXSCHEMENAME      = FluxSchemeName;           // Flux Scheme Name
    Parameters::TIMESCHEMENAME      = TimeSchemeName;           // Time Scheme Name
    Parameters::PHYSICALMODELNAME   = PhysicalModelName;        // Physical model name
    Parameters::SOURCETERMNAME      = SourceTermName;           // Source Term name
    
    Parameters::INITIALFIELD.resize(nbEqs*nbCells);
    for (int i = 0; i < nbEqs*nbCells; i++){
        py::float_ InitialField_py_i = InitialField_py[i];
        Parameters::INITIALFIELD[i] = InitialField_py_i;
    }
}
