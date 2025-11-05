//
//  Parameters.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/07/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "Parameters.hpp"
#include <string>
#include <math.h>

using namespace std;

namespace Parameters {
    int MPI_WRITER             = 0;                    // Root for the MPI
    int NBEQS                  = 1;                    // number of equations
    int NBCELLS                = 100;                  // number of cells
    double DELTAT              = 0.;                   // We make accessible the delta t
    vector<int> NBCELLS_MPI(1);                        // vector with number of cells per processor
    string GEOMETRY            = "1D";                 // geometry, e.g., 1D, 2D, 3D.
    vector<double> MESH(NBCELLS + 1);                      // Vector with the mesh
    double LENGTH              = 1.;                   // Length of the domain
    string INLETTYPE           = "Neumann";            // Type of inlet condition
    string OUTLETTYPE          = "Neumann";            // Type of outlet condition
    vector<double> UINLET(NBEQS);           // solution in the inlet
    vector<double> UOUTLET(NBEQS);          // solution in the outlet
    double A                   = 1.;                   // advection speed of the model
    double CFL                 = 0.45;                 // Courant Friedrichs Lewy number
    double EPS                 = -30.;                  // residual threshold for convergence
    int NBSTEPS                 = 1e9;                 // nb Steps maximum for the time iteration
    int SAVERATE               = 5;                    // save rate
    int PRINTRATE              = 2000;                 // rate to print on screen
    string RESULTDIRECTORY     = "./Results";          // result directory
    string LIMITERNAME         = "ospre";              // Limiter Name
    string RECONSTRUCTORNAME   = "TVD2ndOrder1D";      // Reconstructor Name
    string FLUXSCHEMENAME      = "LaxFriedrich";       // Flux Scheme Name
    string TIMESCHEMENAME      = "ForwardEuler";       // Time Scheme Name
    string SOURCETERMNAME      = "NullSourceTerm";     // Source Term Name
    string PHYSICALMODELNAME   = "AdvectionEq1D";      // Physical Model Name
    string DATAWRITER          = "DataWriter1D";       // Data writer
    double GAMMA               = 5./3.;                // Adiabatic constant
    vector<double> INITIALFIELD(NBEQS*NBCELLS);        // Initial Field
    unsigned int NBFLUIDS      = 1;                    // Number of fluids
    // Options for the two-fluid source term
    double MASSRATIO           = 1.;                   // Mi/Me
    double DEBYELENGTH         = 1.;                   // Normalized Debye lengt
    double IONIZCONST          = 0.;                   // Normalized Ionization Const
    double EPSIONIZ            = 1.;                   // Normalized ionization potential
    double COLLIONS            = 0.;                   // Normalized collisional cross section for ions
    double COLLELECTRONS       = 0.;                   // Normalized collisional cross section for electrons
    double COLLNEGIONS         = 0.;                   // Normalized collisional cross section for negative ions
    vector<double> Mach_in(1);                         // Vector of Mach nb. for which the i-n collision rates has been computed
    vector<double> ListT(1);                           // Vector of temperatures for which the i-n collision rates has been computed
    vector<double> arrayK_in(1);                         // Vector of normalized i-n collision rates (at fixed temperature)
    vector<vector<double>> arrayK_uT(1,vector<double>(1)); // array of normalized i-n momentum collision rates as a function of u and T
    vector<vector<double>> arrayK2_uT(1,vector<double>(1));// array of normalized i-n energy collision rates as a function of u and T
    vector<double> Mach_en(1);                         // Vector of Mach nb. for which the e-n collision rates has been computed
    vector<double> arrayK_en(1);                       // Vector of normalized e-n collision rates (at fixed temperature)
//    double K_iz                = 1.;                 // Ionization rate
    double Nu_iz                = 1.;                  // Ionization frequency
    double PHIIN               = 0.;                   // Normalized boundary for Phi (inlet)
    double PHIOUT              = 0.;                   // Normalized boundary for Phi (outlet)
    vector<double> PHIINITIAL(NBCELLS);                // Vector with the initial value of Phi
    double ATTCONST            = 0.;                   // Attachment constant for three fluids
    double RECCONST            = 0.;                   // Recombination constant for three fluids
    vector<double> MASSES(NBFLUIDS);                   // Normalized masses for multi-fluid models
    double ENTIME              = 0.;                   // Time at which the negative ions are injected
    double PGAS                = 1.;                   // Pressure of the gas in mTorr
    double TGAS                = 300.;                 // Temperature of the gas in K
    double N_REF               = 1e15;                 // Reference density of the plasma
    double T_REF               = 3;                    // Reference electron temperature in eV
    double TIME_REF            = 1e-5;                 // Reference time in seconds
    double L_REF               = 0.03;                 // Reference length in meters
    double PABS                = 1e-6;                 // Normalized absorbed power
    double LY                  = 0.03;                 // length of reactor in y direction in m
    double LZ                  = 0.03;                 // length of reactor in z direction in m
    double OMEGA_RF            = 1403.9856786769978;   // frequency rf normalized
    bool   SHEATHMODEL         = false;                // Flag for sheath mode
    double EPSILON             = 1e-4;                 // Small parameter in Euler Friction
    int STRANG_SPLITTING_STEP  = 0.;                   // Strang-splitting step
    double VOLTAGE             = 200;                  // CCP Voltage
    string TRANSPORTMODEL      = "Maxwellian";         // Transport model for the CCP discharge
    double TWALL               = 300.;                 // Temperature of the wall in K
    vector<double> EXTERNAL_E_FIELD(NBCELLS);          // External electric field
    double INITIAL_TIME        = 0.;                   // Initial time when restarting
    double FINAL_TIME          = 0.;                   // Final time to stop simulation
    double R_S                 = 0.05;                 // Rs for streamer in 1.5D cylindrical
    string GEOMETRY_STREAMER   = "Cylinder";           // Geometry for the streamer





    // Option for the two fluid Isothermal model
    vector<double> SOUNDSPEED(NBFLUIDS);               // vector with the sound speed of the fluids
    vector<double> POLYTROPICINDEX(NBFLUIDS);          // vector with the polytropic indices in the isentropic flow
    vector<double> POLYTROPICCONST(NBFLUIDS);          // vector with the polytropic constants
    vector<double> ETHRUSTERINITIAL(NBCELLS);          // Vector with the constant Electric field for the thruster model
    double CURRENTTHRUSTER     = 0.;                   // current entering the hall thruster
    
    // Constants
    double BOLTZMANN           = 1.380649e-23;         // Boltzmann constant
    double PI                  = atan(1)*4;            // Pi
    double E_CHARGE            = 1.602176634e-19;      // Electron charge

    // Options for Iodine Isothermal
    double GAMMA_WALL_BC  = 0.;         // Recombination at the wall (at x = 0 and x=L)
    double GAMMA_WALL_VOL = 0.;         //  Recombination at the wall (inside volume)
    vector<double> NE_COEFFS(6);        // Coefficients for fitting the temperature
    vector<double> TE_COEFFS(6);        // Coefficients for fitting the temperature

}


