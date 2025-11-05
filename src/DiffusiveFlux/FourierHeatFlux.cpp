//
//  FourierHeatFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FourierHeatFlux.hpp"

void FourierHeatFlux::computeDiffusiveStep(const double dt)
{
    vector<double>& x           = MeshData::getInstance().getData<double>("x");


    // Compute the inner cells
    for (int iCell = 1; iCell < NBCELLS - 1; ++iCell) {
        
        const double e           = 1.602176634e-19;
        //i Cell
        const double rho_i         = (*m_cells)[iCell].uCC[0];
        const double ionMom_z_i    = (*m_cells)[iCell].uCC[1];
        const double TotalEnergy_i = (*m_cells)[iCell].uCC[2];
        const double n_i        = rho_i/m_mi;
        const double x_i        = x[iCell];
        const double dx_i       = (*m_cells)[iCell].dx;
        const double u_x_i      = ionMom_z_i/(rho_i);
        const double p_e_i      = 2./3.*(TotalEnergy_i - 0.5*m_mi*n_i*u_x_i*u_x_i);
        const double T_e_i      = p_e_i/(e*n_i);

        //i-1 Cell
        const double rho_im1    = (*m_cells)[iCell-1].uCC[0];
        const double n_im1      = rho_im1/m_mi;
        const double x_im1      = x[iCell-1];
        //i+1 Cell
        const double rho_ip1    = (*m_cells)[iCell+1].uCC[0];
        const double n_ip1      = rho_ip1/m_mi;
        const double x_ip1      = x[iCell+1];
        
        const double kappa_i      = getConductivity((*m_cells)[iCell].uCC);
        const double kappa_im1    = getConductivity((*m_cells)[iCell - 1].uCC);
        const double kappa_im12   = (kappa_i + kappa_im1)/2;
        const double kappa_ip1    = getConductivity((*m_cells)[iCell + 1].uCC);
        const double kappa_ip12   = (kappa_i + kappa_ip1)/2;
        
        const double factor = 2*dt/(3*n_i*e);
        m_A[iCell][iCell] = 1 + (kappa_ip12/(x_ip1 - x_i) + kappa_im12/(x_i - x_im1))*factor/dx_i;
        m_A[iCell][iCell - 1] = -(kappa_im12/(x_i - x_im1))*factor/dx_i;
        m_A[iCell][iCell + 1] = -(kappa_ip12/(x_ip1 - x_i))*factor/dx_i;
        
        m_B[iCell] = T_e_i;
    }
    // Compute the first cells
    const double e           = 1.602176634e-19;
    //i Cell
    double rho_i         = (*m_cells)[0].uCC[0];
    double ionMom_z_i    = (*m_cells)[0].uCC[1];
    double TotalEnergy_i = (*m_cells)[0].uCC[2];
    double n_i           = rho_i/m_mi;
    double x_i           = x[0];
    double dx_i          = (*m_cells)[0].dx;
    double u_x_i         = ionMom_z_i/(rho_i);
    double p_e_i         = 2./3.*(TotalEnergy_i - 0.5*m_mi*n_i*u_x_i*u_x_i);
    double T_e_i         = p_e_i/(e*n_i);

    //i-1 Cell
    double rho_im1          = m_uInlet[0];
    double ionMom_z_im1     = m_uInlet[1];
    double TotalEnergy_im1  = m_uInlet[2];
    double n_im1            = rho_im1/m_mi;
    double x_im1            = x[0] - dx_i/2.;
    double u_x_im1          = ionMom_z_im1/(rho_im1);
    double p_e_im1          = 2./3.*(TotalEnergy_im1 - 0.5*m_mi*n_im1*u_x_im1*u_x_im1);
    double T_e_im1          = p_e_im1/(e*n_im1);

    //i+1 Cell
    double rho_ip1    = (*m_cells)[1].uCC[0];
    double n_ip1      = rho_ip1/m_mi;
    double x_ip1      = x[1];
    
    double kappa_i      = getConductivity((*m_cells)[0].uCC);
    double kappa_im1    = getConductivity(m_uInlet);
    double kappa_im12   = (kappa_i + kappa_im1)/2;
    double kappa_ip1    = getConductivity((*m_cells)[1].uCC);
    double kappa_ip12   = (kappa_i + kappa_ip1)/2;
    
    double factor = 2*dt/(3*n_i*e);
    m_A[0][0] = 1 + (kappa_ip12/(x_ip1 - x_i) + kappa_im12/(x_i - x_im1))*factor/dx_i;
    m_A[0][1] = -(kappa_ip12/(x_ip1 - x_i))*factor/dx_i;
    m_B[0] = T_e_i + T_e_im1*(kappa_im12/(x_i - x_im1))*factor/dx_i;
    
    // Last cell
    //i Cell
    const int n1 = NBCELLS - 1;
    rho_i         = (*m_cells)[n1].uCC[0];
    ionMom_z_i    = (*m_cells)[n1].uCC[1];
    TotalEnergy_i = (*m_cells)[n1].uCC[2];
    n_i        = rho_i/m_mi;
    x_i        = x[n1];
    dx_i       = (*m_cells)[n1].dx;
    u_x_i      = ionMom_z_i/(rho_i);
    p_e_i      = 2./3.*(TotalEnergy_i - 0.5*m_mi*n_i*u_x_i*u_x_i);
    T_e_i      = p_e_i/(e*n_i);

    //i-1 Cell
    rho_im1    = (*m_cells)[n1-1].uCC[0];
    n_im1      = rho_im1/m_mi;
    x_im1      = x[n1-1];
    //i+1 Cell
    rho_ip1                 = m_uOutlet[0];
    double ionMom_z_ip1     = m_uOutlet[1];
    double TotalEnergy_ip1  = m_uOutlet[2];
    n_ip1                   = rho_ip1/m_mi;
    x_ip1                   = x[n1] + dx_i/2.;
    double u_x_ip1          = ionMom_z_ip1/(rho_ip1);
    double p_e_ip1          = 2./3.*(TotalEnergy_ip1 - 0.5*m_mi*n_ip1*u_x_ip1*u_x_ip1);
    double T_e_ip1          = p_e_ip1/(e*n_ip1);
    
    kappa_i      = getConductivity((*m_cells)[n1].uCC);
    kappa_im1    = getConductivity((*m_cells)[n1 - 1].uCC);
    kappa_im12   = (kappa_i + kappa_im1)/2;
    kappa_ip1    = getConductivity(m_uOutlet);
    kappa_ip12   = (kappa_i + kappa_ip1)/2;
    
    factor = 2*dt/(3*n_i*e);
    m_A[n1][n1] = 1 + (kappa_ip12/(x_ip1 - x_i) + kappa_im12/(x_i - x_im1))*factor/dx_i;
    m_A[n1][n1 - 1] = -(kappa_im12/(x_i - x_im1))*factor/dx_i;
    m_B[n1] = T_e_i + T_e_ip1*(kappa_ip12/(x_ip1 - x_i))*factor/dx_i;

    // Compute the variable Phi (that is in m_x)
    m_linearSolver->solveLinearSystem(m_A, m_x, m_B);

    for (unsigned int iCell = 0; iCell < NBCELLS; ++iCell){
        const double rho_i         = (*m_cells)[iCell].uCC[0];
        const double ionMom_z_i    = (*m_cells)[iCell].uCC[1];
        const double TotalEnergy_i = (*m_cells)[iCell].uCC[2];
        const double n_i        = rho_i/m_mi;
        const double x_i        = x[iCell];
        const double dx_i       = (*m_cells)[iCell].dx;
        const double u_x_i      = ionMom_z_i/(rho_i);
        const double p_e_i      = 2./3.*(TotalEnergy_i - 0.5*m_mi*n_i*u_x_i*u_x_i);
        const double T_e_i      = p_e_i/(e*n_i);

        (*m_cells)[iCell].uCC[2] = 3./2.*m_x[iCell]*n_i*e + 0.5*rho_i*u_x_i*u_x_i;
    }
}
