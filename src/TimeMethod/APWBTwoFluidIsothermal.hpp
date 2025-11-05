//
//  APWBTwoFluidIsothermal.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 21/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef APWBTwoFluidIsothermal_hpp
#define APWBTwoFluidIsothermal_hpp


#include <stdio.h>
#include <vector>
#include <iostream>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"
#include "../SpaceMethod/SpaceMethod.hpp"
#include "../FluxScheme/FluxScheme.hpp"
#include "../SourceTerm/SourceTerm.hpp"
#include "../SpaceReconstructor/SpaceReconstructor.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"
#include "../PhysicalModel/PhysicalModel.hpp"
#include "../LinearSolver/LinearSolver.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../MeshData/MeshData.hpp"
#include <cmath>

using namespace std;
using namespace Parameters;


class APWBTwoFluidIsothermal : public TimeMethod
{
public:
    APWBTwoFluidIsothermal(string name) : TimeMethod(name), u_1(NBEQS*NBCELLS), m_A(NBCELLS,vector<double>(NBCELLS,0.)), m_B(NBCELLS, 0.), m_gradPhi(NBCELLS,0.), m_inverseDensity(NBCELLS,1.), m_density_nP1Minus(NBCELLS, 0.), m_momentum_nP1Minus(NBCELLS, 0.),m_density_nP1(NBCELLS, 0.), m_momentum_nP1(NBCELLS, 0.), m_velocity(NBCELLS, 0.), m_localCell_L(2,0.), m_localCell_R(2,0.){
        cout<<"Entering "<< this->getName() <<"\n";
        setup();
        cout<<"After setup\n";
        
        // Initialize the matrix
        // Set up initial field of Phi
        m_x = Parameters::PHIINITIAL;
        
        // Set up the matrix
        m_A[0][0] = -2;
        m_A[0][1] = 1;
        unsigned int n1 = NBCELLS - 1;
        // Inner cells
        for (unsigned int i = 1; i < n1; i++) {
            m_A[i][i]     = -2;
            m_A[i][i - 1] = 1;
            m_A[i][i + 1] = 1;
        }
        // Last cell
        m_A[n1][n1]   = -2;
        m_A[n1][n1-1] = 1;
    }
    ~APWBTwoFluidIsothermal() {}
    void setup();
    void spatialDiscretization(){setBoundaries(); computeIonFlux();}
    void computeElectronDensity();
    void computeElectricPotential();
    void computeElectronVelocity();
    void computeEulerianStep();
    double computeIonizationConstant();
    void computeSources();
    double velocityFlux(const CellDataRef uL, const CellDataRef uR);
    double velocityFlux_Euler(const CellDataRef uL, const CellDataRef uR);
    double pressureFlux(const CellDataRef uL, const CellDataRef uR);
    double densityNumericalViscosity(const double niP1,const double ni, const double niM1);
    void setBoundaries();
    void computeIonFlux();
    void setDt(){
        vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
        const double Dx = cells[0].dx;
        double n_e, u_e, lambda_max, omega_max = 0;
        const double c_e      = sqrt(MASSRATIO);
        m_Dt = CFL/sqrt(cells[0].uCC[0])*DEBYELENGTH*sqrt(MASSRATIO); // Initialization of the Dt
        for (int iCell = 0; iCell < NBCELLS; ++iCell) {
            n_e = cells[iCell].uCC[0];
            u_e = cells[iCell].uCC[1];
            lambda_max = max(abs(u_e - c_e), abs(u_e + c_e));
            omega_max  = sqrt(n_e)*DEBYELENGTH*sqrt(MASSRATIO);
            
            m_Dt = min(m_Dt, min(CFL*Dx/lambda_max, CFL/omega_max));
            
        }
        //cout<<"Dt = "<< m_Dt <<"\n";
    }
    double getDtOvDx(){
        const double Dx = getDx();
        
        return m_Dt/Dx;
    }
    
    double getDx() {return m_Dx;}
    void setDx(const double Dx) {m_Dx = Dx;}
    virtual void takeStep(double dt);
    void takeStep();
    double WBmassSourceL(const double u_12, const double c_12);
    double WBmassSourceR(const double u_12, const double c_12);
    double WBmomSourceL(const double u_12, const double c_12);
    double WBmomSourceR(const double u_12, const double c_12);

    double signfunction(const double x){
//        double a = 0;
//        if (x > 0){a = 1;}
//        else if (x < 0){ a = -1;}
//        else {a = 0;}
//        return a;
        return tanh(x/(0.1*SOUNDSPEED[1]));
        
    }
    
protected:
    vector<double> u_1;
    double m_dtOvdx;
    // Needed for the space discretization
    unique_ptr<FluxScheme> m_flux;
    unique_ptr<SourceTerm> m_source;
    unique_ptr<SpaceReconstructor> m_reconstructor;
    unique_ptr<BoundaryCondition> m_InletBC;
    unique_ptr<BoundaryCondition> m_OutletBC;
    unique_ptr<PhysicalModel> m_pm;
    unique_ptr<LinearSolver> m_linearSolver;
    Matrix m_A;         // Matrix of the linear solver
    vector<double> m_x; // Solution of the linear solver
    vector<double> m_B; // rhs of the linear solver
    vector<double> m_gradPhi; // gradient of the electric potential
    
    
    CellDataRef m_uInlet;
    CellDataRef m_uOutlet;
    vector<double> m_Fip12;
    vector<double> m_Fim12;
    vector<double> m_inverseDensity;
    vector<double> m_density_nP1Minus;
    vector<double> m_momentum_nP1Minus;
    vector<double> m_density_nP1;
    vector<double> m_momentum_nP1;
    vector<double> m_velocity;
    
    vector<double> m_localCell_L;
    vector<double> m_localCell_R;
    
    double m_Dx;
    double m_Dt;
};

#endif /* APWBTwoFluidIsothermal_hpp */
