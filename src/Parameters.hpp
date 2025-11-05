//
//  Parameters.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/07/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Parameters_hpp
#define Parameters_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


namespace py = pybind11;
using namespace std;

namespace Parameters {

    /// global variables (simulation invariants)
    extern int MPI_WRITER;            // root for the MPI
    extern int NBEQS;               // number of equations
    extern int NBCELLS;             // number of cells
    extern double DELTAT;           // We make accessible the delta t
    extern vector<int> NBCELLS_MPI; // number of cells
    extern string GEOMETRY;         // geometry, e.g., 1D, 2D, 3D.
    extern vector<double> MESH;     // Vector with the mesh
    extern double LENGTH;           // Length of the domain
    extern string INLETTYPE;        // Type of inlet condition
    extern string OUTLETTYPE;       // Type of outlet condition
    extern vector<double> UINLET;   // solution in the inlet
    extern vector<double> UOUTLET;  // solution in the outlet
    extern double A;                // advection speed of the model
    extern double CFL;              // Courant Friedrichs Lewy number
    extern double EPS;              // residual threshold for convergence
    extern int NBSTEPS;             // nb Steps maximum for the time iteration
    extern int SAVERATE;            // rate to save the result
    extern int PRINTRATE;           // rate to print on screen
    extern string RESULTDIRECTORY;  // result directory
    extern string LIMITERNAME;      // Limiter name
    extern string RECONSTRUCTORNAME;// Reconstructor name
    extern string FLUXSCHEMENAME;   // Flux Scheme Name
    extern string TIMESCHEMENAME;   // Time Scheme Name
    extern string SOURCETERMNAME;   // Source Term Name
    extern vector<double> INITIALFIELD; // Initial field
    extern string PHYSICALMODELNAME;    // Physical Model Name
    extern string DATAWRITER;           // Data Writer
    extern double GAMMA;                // Adiabatic ratio
    extern unsigned int NBFLUIDS;       // Number of fluids
    // Options for the two fluids source term
    extern double MASSRATIO;        // Mi/Me
    extern double DEBYELENGTH;      // Normalized Debye lengt
    extern double IONIZCONST;       // Normalized Ionization Const
    extern double EPSIONIZ;         // Normalized ionization potential
    extern double COLLIONS;         // Normalized collisional cross section for ions
    extern double COLLELECTRONS;    // Normalized collisional cross section for electrons
    extern double COLLNEGIONS;      // Normalized collisional cross section for negative ions
    extern vector<double> Mach_in; // Vector of Mach nb. for which the i-n coll. rates has been computed
    extern vector<double> arrayK_in; // Vector of normalized i-n collision rates (at fixed temperature)
    extern vector<double> ListT; // Vector of temperature for which the i-n coll. rates has been computed
    extern vector<vector<double>> arrayK_uT; // array of normalized i-n momentum collision rates as a function of u and T
    extern vector<vector<double>> arrayK2_uT; // array of normalized i-n energy collision rates as a function of u and T
    extern vector<double> Mach_en; // Vector of Mach nb. for which the e-n coll. rates has been computed
    extern vector<double> arrayK_en; // Vector of normalized e-n collision rates (at fixed temperature)
//    extern double K_iz;             // Ionization rate
    extern double Nu_iz;            // Ionization frequency
    extern double PHIIN;            // Normalized boundary for Phi (inlet)
    extern double PHIOUT;           // Normalized boundary for Phi (outlet)
    extern vector<double> PHIINITIAL; // Vector with the initial value of Phi
    extern vector<double> ETHRUSTERINITIAL; // Vector with the constant Electric field for the thruster model
    extern double ATTCONST;                 // Attachment constant for three fluids
    extern double RECCONST ;                // Recombination constant for three fluids
    extern vector<double> MASSES;           // Normalized masses for multi-fluid models
    extern double ENTIME;                   // Time at which the negative ions are injected
    extern double PGAS;                     // Pressure of the gas in mTorr
    extern double TGAS;                     // Temperature of the gas in K
    extern double N_REF;                    // Reference density of the plasma
    extern double T_REF;                    // Reference electron temperature in eV
    extern double TIME_REF;                 // Reference time in seconds
    extern double L_REF;                    // Reference length in meters
    extern double PABS;                     // Normalized absorbed power
    extern double LY;                       // length of reactor in y direction in m
    extern double LZ;                       // length of reactor in z direction in m
    extern double OMEGA_RF;                 // length of reactor in z direction in m
    extern bool   SHEATHMODEL;              // Flag for sheath model
    extern double EPSILON;                  // Small parameter in Euler Friction
    extern int STRANG_SPLITTING_STEP;       // Strang-splitting step
    extern double VOLTAGE;                  // CCP Voltage
    extern string TRANSPORTMODEL;          // Transport model for the CCP discharge
    extern double TWALL;                 // Temperature of the wall in K
    extern vector<double> EXTERNAL_E_FIELD;        // External electric field
    extern double INITIAL_TIME;            // Initial time when restarting
    extern double FINAL_TIME;              // Final time to stop simulation
    extern double R_S;                     // Final time to stop simulation
    extern string GEOMETRY_STREAMER;                // Streamer geometry

   

    
    // Option for the two fluid Isothermal model
    extern vector<double> SOUNDSPEED; // vector with the sound speed of the fluids
    extern vector<double> POLYTROPICINDEX;          // vector with the polytropic indices in the isentropic flow
    extern vector<double> POLYTROPICCONST;          // vector with the polytropic constants
    extern double CURRENTTHRUSTER;                  // current entering the hall thruster
    
    // Constants
    extern double BOLTZMANN;                        // Boltzmann constant
    extern double PI;                               // Pi
    extern double E_CHARGE;                         // Electron charge

    // Options for Iodine Isothermal
    extern double GAMMA_WALL_BC;                    // Recombination at the wall (at x = 0 and x=L)
    extern double GAMMA_WALL_VOL;                   //  Recombination at the wall (inside volume)
    extern vector<double> NE_COEFFS;                // Coefficients for fitting the temperature
    extern vector<double> TE_COEFFS;                // Coefficients for fitting the temperature
    
}

#endif /* Parameters_hpp */
