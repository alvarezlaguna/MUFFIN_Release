//
//  SingleIonMixture.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SingleIonMixture_hpp
#define SingleIonMixture_hpp

#include <stdio.h>
#include "Mixture.hpp"
#include <mpi.h>

using namespace Parameters;

class SingleIonMixture : public Mixture
{
public:
    SingleIonMixture(string name) : Mixture(name), m_K_rates(), m_Full_Losses(), m_Omega11_en(), m_Omega12_en(), m_Omega13_en(), m_Omega14_en(), m_Omega22_ee(), m_Omega23_ee(), m_Omega24_ee(), m_K_ionization(), m_IonOmegaM_rates(), m_IonOmegaE_rates(), m_LossesFourthDelta(), m_LossesDelta(), m_IzDelta(), m_Temperature(), m_Delta(), m_IonMach() {}
    ~SingleIonMixture(){}
    
    void setup();
    void unsetup();
    
    void setIonMass(double mIon) { m_IonMass = mIon;}
    void setElecMass(double mElec) { m_ElecMass = mElec;}
    void setGasDensity(double nGas) {m_n_gas = nGas;}
    void setGasTemperature(double TGas) {m_T_gas = TGas;}
    double getIonMass() {return m_IonMass;}
    double getElecMass() {return m_ElecMass;}
    double getGasTemperature() {return m_T_gas;}
    double getGasDensity() {return m_n_gas;}
    double computeIonizationRate(const double& Te);
    double computeIonizationRate(const double& Te, const double& Deltae);
    double computeElecInelEnergyLosses(const double& Te);
    double computeElecInelEnergyLosses(const double& Te, const double& Deltae);
    double computeElecInelFourthMomLosses(const double& Te, const double& Deltae);
    double computeOmega11_en(const double& Te);
    double computeOmega12_en(const double& Te);
    double computeOmega13_en(const double& Te);
    double computeOmega14_en(const double& Te);
    double computeOmega22_ee(const double& Te);
    double computeOmega23_ee(const double& Te);
    double computeOmega24_ee(const double& Te);
    
    double computeElecRate(const double Te,const vector<double>& coeffs){
        
//        if(Te < 2){
//            int my_rank;
//            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//            if(my_rank == MPI_WRITER){
//                py::print("\n**************************************\n");
//                py::print("\nWARNING: Fit for the rate might be wrong as Te < 2 eV");
//                py::print("\n**************************************\n\n\n");}
//
//        }
        
        return coeffs[0]*pow(Te,(-coeffs[5]))*exp(-coeffs[1]/Te + coeffs[2]/pow(Te, 2) - coeffs[3]/pow(Te, 3)  + coeffs[4]/pow(Te, 4));
    }
    
    double computeIonRate(const vector<double>& Coeffs, const double T){

        const double a = Coeffs[0];
        const double b = Coeffs[1];
        const double c = Coeffs[2];
        const double d = Coeffs[3];

        
        const double t = (1/T);
        
        return a*pow(t,d)*exp(-b*t+ c*pow(t,2));
    }
    
    double linearInterpolation(const double x0, const double x1, const double y0, const double y1, const double x){
        return (y0*(x1 - x) + y1*(x - x0))/(x1 - x0);
    }
    const vector<double> &getOmegaMCoeffs(const int i_Coll, const int i_Mach){
        return m_IonOmegaM_rates[i_Coll*m_NumberMachEntries + i_Mach];
    }
    const vector<double> &getOmegaECoeffs(const int i_Coll, const int i_Mach){
        return m_IonOmegaE_rates[i_Coll*m_NumberMachEntries + i_Mach];
    }
    const double getMach(const int i_Coll, const int i_Mach){
        return m_IonMach[i_Coll*m_NumberMachEntries + i_Mach];
    }
    const double getDelta(const int i_Delta){
        return m_Delta[i_Delta];
    }
    const double getTemperature(const int i_Temp){
        return m_Temperature[i_Temp];
    }
    double computeIonOmegaM(const double& T, const double& Mach);
    double computeIonOmegaE(const double& T, const double& Mach);
    
    
protected:
    vector<vector<double>> m_K_rates;
    vector<double> m_Full_Losses;
    vector<double> m_Omega11_en;
    vector<double> m_Omega12_en;
    vector<double> m_Omega13_en;
    vector<double> m_Omega14_en;
    vector<double> m_Omega22_ee;
    vector<double> m_Omega23_ee;
    vector<double> m_Omega24_ee;
    vector< vector<double> > m_IonOmegaM_rates;
    vector< vector<double> > m_IonOmegaE_rates;
    vector< vector<double> > m_LossesFourthDelta;
    vector< vector<double> > m_LossesDelta;
    vector< vector<double> > m_IzDelta;
    vector< double > m_Temperature;
    vector< double > m_Delta;
    vector<double> m_IonMach;
    vector<double> m_K_ionization;
    double m_IonMass;
    double m_ElecMass;
    double m_n_gas;
    double m_T_gas;
    int m_NumberMachEntries;
    int m_NumberIonColls;
    double m_DeltaMach;

    double m_charge;

};

#endif /* SingleIonMixture_hpp */
