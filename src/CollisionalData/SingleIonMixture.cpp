//
//  SingleIonMixture.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "SingleIonMixture.hpp"
#include "CollisionalData.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void SingleIonMixture::setup(){
    // Set-up the data from the cross sections of the processes
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(my_rank == MPI_WRITER)
        py::print("Setting-up Mixture class:\t ",this->getName(),"\n");
    
    // Get the collisional cross-section data
    CollisionalData& cd = CollisionalData::getInstance();
    auto& mixData  = cd.getData<py::handle>("mixtureData");
    auto& elecData = cd.getData<py::handle>("electronData");
    auto& ionData  = cd.getData<py::handle>("ionData");
    
    // Loop over the collisional processes of electrons
    int counter_collision = 0;
    for (auto collision : elecData){
        if(my_rank == MPI_WRITER)
            py::print("\t Collision Type  : ", *collision.attr("name"));
        // Computing the integration over the energy space
        if(my_rank == MPI_WRITER)
            py::print("\t \t Computing the rate");
        *collision.attr("computeRate")();
        // Computing the integration over the coefficients
        if(my_rank == MPI_WRITER)
            py::print("\t \t Computing coefficients for fit");
        *collision.attr("fitRate")();
        string nameCol = py::str(*collision.attr("name"));
        nameCol = nameCol.c_str(); //conversion needed
        // I save the ionization in a separate vector to avoid loosing time with find
        if(nameCol.find("Ioniz") != string::npos){
            for (auto Coeffs : *collision.attr("Coeffs")){
                double Coefficient = Coeffs.cast<float>();
                m_K_ionization.push_back(Coefficient);
            }
        }
        // Save the Ionization with Deltae
        if(nameCol.find("Ioniz") != string::npos){
            vector<double> Coefficients(1000);
            int i_Coeffs;
            int i_Delta = 0;
            *collision.attr("computeRateDeltaZero")();
            for (auto array : *collision.attr("KDeltaZero")){
                i_Coeffs = 0;
                for(auto Coeffs : array){
                    double Coefficient = Coeffs.cast<float>();
                    Coefficients[i_Coeffs] = Coefficient;
                    i_Coeffs++;
                }
                m_IzDelta.push_back(Coefficients);
                i_Delta++;
            }
        }
        // I save all the coefficients into a map with the potential in the end.
        vector<double> Coefficients(7);
        int i_coeff = 0;
        for (auto Coeffs : *collision.attr("Coeffs")){
            Coefficients[i_coeff] = Coeffs.cast<float>();
            ++i_coeff;
        }
        
        // We put the excitation/ionization potential in the last position
        auto potential_py = *collision.attr("potential");
        double potential  = potential_py.cast<float>();
        Coefficients[6] = potential;
        m_K_rates.push_back(Coefficients);
        Coefficients.clear();
        
        // Test
        *collision.attr("computeKInelasticApprox")(3., 0.1, 0.);
        
        // Read the charge from scipy
        m_charge = 1.602176634e-19;
        ++counter_collision;
        
    }

    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the full losses");
    *mixData.attr("computeFullLosses")();
    *mixData.attr("fitFullLosses")();
    // Compute Full Losses Fourth Delta
    *mixData.attr("computeFullLossesFourthDelta")();
    int i_Coeffs;
    int i_Delta = 0;
    vector<double> Coefficients(1000);
    for (auto array : *mixData.attr("LossesFourthDelta")){
        i_Coeffs = 0;
        for(auto Coeffs : array){
            double Coefficient = Coeffs.cast<float>();
            Coefficients[i_Coeffs] = Coefficient;
            i_Coeffs++;
        }
        m_LossesFourthDelta.push_back(Coefficients);
        i_Delta++;
    }
    // Compute Full Losses Delta
    *mixData.attr("computeFullLossesDelta")();
    i_Delta = 0;
    for (auto array : *mixData.attr("LossesDelta")){
        i_Coeffs = 0;
        for(auto Coeffs : array){
            double Coefficient = Coeffs.cast<float>();
            Coefficients[i_Coeffs] = Coefficient;
            i_Coeffs++;
        }
        m_LossesDelta.push_back(Coefficients);
        i_Delta++;
    }
    for (auto delta_py : *mixData.attr("Delta")){
        double delta = delta_py.cast<float>();
        m_Delta.push_back(delta);
    }
    for (auto temp_py : *mixData.attr("Temperature")){
        double temp = temp_py.cast<float>();
        m_Temperature.push_back(temp);
    }
    
    auto FullLossesCoeffs = *mixData.attr("CoeffsLosses");
    for(auto Coeff : FullLossesCoeffs){
        double Full_Losses_Coeff = Coeff.cast<float>();
        m_Full_Losses.push_back(Full_Losses_Coeff);
    }
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega11_en");
    *mixData.attr("computeOmega11_en")();
    *mixData.attr("fitOmega11_en")();
    auto Omega11_enCoeffs = *mixData.attr("CoeffsOmega11_en");
    for(auto Coeff : Omega11_enCoeffs){
        double Omega11_en = Coeff.cast<float>();
        m_Omega11_en.push_back(Omega11_en);
    }
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega12_en");
    *mixData.attr("computeOmega12_en")();
    *mixData.attr("fitOmega12_en")();
    auto Omega12_enCoeffs = *mixData.attr("CoeffsOmega12_en");
    for(auto Coeff : Omega12_enCoeffs){
        double Omega12_en = Coeff.cast<float>();
        m_Omega12_en.push_back(Omega12_en);
    }
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega13_en");
    *mixData.attr("computeOmega13_en")();
    *mixData.attr("fitOmega13_en")();
    auto Omega13_enCoeffs = *mixData.attr("CoeffsOmega13_en");
    for(auto Coeff : Omega13_enCoeffs){
        double Omega13_en = Coeff.cast<float>();
        m_Omega13_en.push_back(Omega13_en);
    }
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega14_en");
    *mixData.attr("computeOmega14_en")();
    *mixData.attr("fitOmega14_en")();
    auto Omega14_enCoeffs = *mixData.attr("CoeffsOmega14_en");
    for(auto Coeff : Omega14_enCoeffs){
        double Omega14_en = Coeff.cast<float>();
        m_Omega14_en.push_back(Omega14_en);
    }
    
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega22_ee");
    *mixData.attr("computeOmega22_ee")();
    *mixData.attr("fitOmega22_ee")();
    auto Omega22_eeCoeffs = *mixData.attr("CoeffsOmega22_ee");
    for(auto Coeff : Omega22_eeCoeffs){
        double Omega22_ee = Coeff.cast<float>();
        m_Omega22_ee.push_back(Omega22_ee);
    }
    
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega23_ee");
    *mixData.attr("computeOmega23_ee")();
    *mixData.attr("fitOmega23_ee")();
    auto Omega23_eeCoeffs = *mixData.attr("CoeffsOmega23_ee");
    for(auto Coeff : Omega23_eeCoeffs){
        double Omega23_ee = Coeff.cast<float>();
        m_Omega23_ee.push_back(Omega23_ee);
    }
    
    if(my_rank == MPI_WRITER)
        py::print("\t Computing Fit for the Omega24_ee");
    *mixData.attr("computeOmega24_ee")();
    *mixData.attr("fitOmega24_ee")();
    auto Omega24_eeCoeffs = *mixData.attr("CoeffsOmega24_ee");
    for(auto Coeff : Omega24_eeCoeffs){
        double Omega24_ee = Coeff.cast<float>();
        m_Omega24_ee.push_back(Omega24_ee);
    }
    
    
    if(my_rank == MPI_WRITER)
        py::print("\n");
    // Loop over the collisional processes of ions
    int i_ionColl = 0;
    for (auto collision : ionData){
        if(my_rank == MPI_WRITER)
            py::print("\t Collision Type  : ", *collision.attr("name"));
        // Computing the integration over the energy space
        if(my_rank == MPI_WRITER)
            py::print("\t \t Computing the Momentum Rate");
        *collision.attr("computeOmegaM")();
        *collision.attr("fitOmegaM")();
        if(my_rank == MPI_WRITER)
            py::print("\t \t Computing the Energy Rate");
        *collision.attr("computeOmegaE")();
        *collision.attr("fitOmegaE")();
        
        string nameCol = py::str(*collision.attr("name"));
        nameCol = nameCol.c_str(); //conversion needed
        for (auto Mach : *collision.attr("Mach")){
            double Mach_fl = Mach.cast<float>();
            m_IonMach.push_back(Mach_fl);
        }
        int i_Mach = 0;
        vector<double> Coefficients(4);
        int i_Coeffs;
        for (auto array : *collision.attr("CoeffsOmegaM")){
            i_Coeffs = 0;
            for(auto Coeffs : array){
                double Coefficient = Coeffs.cast<float>();
                Coefficients[i_Coeffs] = Coefficient;
                i_Coeffs++;
            }
            m_IonOmegaM_rates.push_back(Coefficients);
            i_Mach++;
        }
        i_Mach = 0;
        for (auto array : *collision.attr("CoeffsOmegaE")){
            i_Coeffs = 0;
            for(auto Coeffs : array){
                double Coefficient = Coeffs.cast<float>();
                Coefficients[i_Coeffs] = Coefficient;
                i_Coeffs++;
            }
            m_IonOmegaE_rates.push_back(Coefficients);
            i_Mach++;

        }
        m_DeltaMach = m_IonMach[1] - m_IonMach[0];
        cout<<"m_Delta_Mach = "<<m_DeltaMach<<"\n";
        m_NumberMachEntries = i_Mach;
        cout<<"\n**************************************\n\n\n";
        cout<<"\n WARNING : We assume that the Mach numbers \n";
        cout<<"   are equally spaced                        \n";
        cout<<"\n**************************************\n\n\n\n";

        
        auto ionMass = *collision.attr("mass");
        
        setIonMass(ionMass.cast<float>());
        setElecMass(9.1093837015e-31);
        ++i_ionColl;
    }
    m_NumberIonColls = i_ionColl;
    
    auto& p_gas  = cd.getData<py::handle>("p_gas");
    auto& T_gas  = cd.getData<py::handle>("T_gas");
    double n_gas = 0.13332237*p_gas.cast<float>()/(E_CHARGE*T_gas.cast<float>());
    setGasDensity(n_gas);
    setGasTemperature(T_gas.cast<float>());
    
    
    //Example to identify one collision type
    //            string nameCol = py::str(*collision.attr("name"));
    //            if(std::strncmp(nameCol.c_str(), "Ar-e (14.00 eV)", 15) == 0){
}

void SingleIonMixture::unsetup(){

}

double SingleIonMixture::computeIonizationRate(const double& Te){
    
    if(m_K_ionization.empty()){
        cout<<"\n**************************************\n\n";
        cout<<"\nERROR: No Ionization in the mixture\n";
        cout<<"\n**************************************\n\n\n\n";
        exit(0);
    }

    return computeElecRate(Te, m_K_ionization);
    
}

double SingleIonMixture::computeIonizationRate(const double& Te, const double& Deltae){
    
    if(m_K_ionization.empty()){
        cout<<"\n**************************************\n\n";
        cout<<"\nERROR: No Ionization in the mixture\n";
        cout<<"\n**************************************\n\n\n\n";
        exit(0);
    }
//    // Get the collisional cross-section data
//    CollisionalData& cd = CollisionalData::getInstance();
//    auto& mixData  = cd.getData<py::handle>("mixtureData");
//    auto& elecData = cd.getData<py::handle>("electronData");
//    double freq = 0.;
//    for (auto collision : elecData){
//
//        string nameCol = py::str(*collision.attr("name"));
//        nameCol = nameCol.c_str(); //conversion needed
//        // I save the ionization in a separate vector to avoid loosing time with find
//        if(nameCol.find("Ioniz") != string::npos){
//            auto freq_py = *collision.attr("computeKInelasticApprox")(Te, Deltae, 0.);
//            freq =freq_py.cast<float>();
//        }
//    }
//    cout<<"Old = "<<freq<<"\n";
    
    // TODO: WE RESTRICT TO DELTA\in(-1,1) and Te in (1,20)
    double freq = 0.;
    double DeltaEff = Deltae;
    double TeEff  = Te;
    if(Deltae < -1.){
        DeltaEff = -1.;
    }
    if(Deltae > 1.){
        DeltaEff = 1.;
    }
    if(Te < 1.){
        TeEff = 1.;
    }
    if(TeEff > 20.){
        TeEff = 20.;
    }
    double Losses = 0.;
    double Delta_Delta = m_Delta[1] - m_Delta[0];
    double Delta_Temp  = m_Temperature[1] - m_Temperature[0];
    int lowerDeltaIdx = (int) (std::floor((DeltaEff - m_Delta[0])/Delta_Delta));
    int upperDeltaIdx = lowerDeltaIdx + 1;
    int lowerTempIdx = (int) (std::floor((TeEff-m_Temperature[0])/Delta_Temp));
    int upperTempIdx = lowerTempIdx + 1;

    double lowerDelta_lowerTemp = m_IzDelta[lowerDeltaIdx][ lowerTempIdx];
    double lowerDelta_upperTemp = m_IzDelta[lowerDeltaIdx][ upperTempIdx];
    double upperDelta_lowerTemp = m_IzDelta[upperDeltaIdx][ lowerTempIdx];
    double upperDelta_upperTemp = m_IzDelta[upperDeltaIdx][ upperTempIdx];
    
    double lower_Delta = linearInterpolation(m_Temperature[lowerTempIdx],m_Temperature[upperTempIdx] , lowerDelta_lowerTemp, lowerDelta_upperTemp, TeEff);
    double upper_Delta = linearInterpolation(m_Temperature[lowerTempIdx],m_Temperature[upperTempIdx] , upperDelta_lowerTemp, upperDelta_upperTemp, TeEff);
    freq = linearInterpolation(m_Delta[lowerDeltaIdx],m_Delta[upperDeltaIdx] , lower_Delta, upper_Delta, DeltaEff);

    return freq;
}

double SingleIonMixture::computeElecInelEnergyLosses(const double& Te){

    // Very slow method that requires a loop over all the processes
//    double losses = 0.;
//    for (auto collision : m_K_rates){
//        double rate          = computeElecRate(Te, collision);
//        double exc_potential = (collision).back(); //Last value of the vector
//
//        losses += (-1)*m_charge*rate*exc_potential;
//
//    }
    
    double losses = (-1)*m_charge*computeElecRate(Te, m_Full_Losses);
    
    return losses;
}

double SingleIonMixture::computeElecInelEnergyLosses(const double& Te, const double& Deltae){
    
//    double losses = 0.;
//
//    CollisionalData& cd = CollisionalData::getInstance();
//    auto& mixData  = cd.getData<py::handle>("mixtureData");
//    auto& elecData = cd.getData<py::handle>("electronData");
//    for (auto collision : elecData){
//
//        string nameCol = py::str(*collision.attr("name"));
//        nameCol = nameCol.c_str(); //conversion needed
//        auto freq_py = *collision.attr("computeKInelasticApprox")(Te, Deltae, 0.);
//        double freq =freq_py.cast<float>();
//        auto potential_py = *collision.attr("potential");
//        double potential  = potential_py.cast<float>();
//        losses += (-1)*m_charge*freq*potential;
//    }
//    cout<<"Old =" <<losses<<"\n";
    
    // TODO: WE RESTRICT TO DELTA\in(-1,1) and Te in (1,20)
    double DeltaEff = Deltae;
    double TeEff  = Te;
    if(Deltae < -1.){
        DeltaEff = -1.;
    }
    if(Deltae > 1.){
        DeltaEff = 1.;
    }
    if(Te < 1.){
        TeEff = 1.;
    }
    if(TeEff > 20.){
        TeEff = 20.;
    }
    double Losses = 0.;
    double Delta_Delta = m_Delta[1] - m_Delta[0];
    double Delta_Temp  = m_Temperature[1] - m_Temperature[0];
    int lowerDeltaIdx = (int) (std::floor((DeltaEff - m_Delta[0])/Delta_Delta));
    int upperDeltaIdx = lowerDeltaIdx + 1;
    int lowerTempIdx = (int) (std::floor((TeEff-m_Temperature[0])/Delta_Temp));
    int upperTempIdx = lowerTempIdx + 1;

    double lowerDelta_lowerTemp = m_LossesDelta[lowerDeltaIdx][ lowerTempIdx];
    double lowerDelta_upperTemp = m_LossesDelta[lowerDeltaIdx][ upperTempIdx];
    double upperDelta_lowerTemp = m_LossesDelta[upperDeltaIdx][ lowerTempIdx];
    double upperDelta_upperTemp = m_LossesDelta[upperDeltaIdx][ upperTempIdx];
    
    double lower_Delta = linearInterpolation(m_Temperature[lowerTempIdx],m_Temperature[upperTempIdx] , lowerDelta_lowerTemp, lowerDelta_upperTemp, TeEff);
    double upper_Delta = linearInterpolation(m_Temperature[lowerTempIdx], m_Temperature[upperTempIdx] , upperDelta_lowerTemp, upperDelta_upperTemp, TeEff);
    Losses = linearInterpolation(m_Delta[lowerDeltaIdx],m_Delta[upperDeltaIdx] , lower_Delta, upper_Delta, DeltaEff);

    return -(1.)*Losses*m_charge;
}

double SingleIonMixture::computeElecInelFourthMomLosses(const double& Te, const double& Deltae){
    
//    double losses = 0.;
//    CollisionalData& cd = CollisionalData::getInstance();
//    auto& mixData  = cd.getData<py::handle>("mixtureData");
//    auto& elecData = cd.getData<py::handle>("electronData");
//    for (auto collision : elecData){
//
//        string nameCol = py::str(*collision.attr("name"));
//        nameCol = nameCol.c_str(); //conversion needed
//        cout<<"nameCol = "<<nameCol<<"\n";
//        auto freq_py = (*collision.attr("computeKInelasticApprox")(Te, Deltae, 0));
//        double freq  = freq_py.cast<float>();
//
//        auto freq2_py = (*collision.attr("computeKInelasticApprox")(Te, Deltae, 1.));
//        double freq2 = freq2_py.cast<float>();
//        cout<<"freq inside code = "<<freq<<"\n";
//        cout<<"freq2 inside code = "<<freq2<<"\n";
//        auto potential_py = *collision.attr("potential");
//        double potential  = potential_py.cast<float>();
//        losses += -freq*potential*potential/(Te*Te) -freq2*potential/Te;
//    }
//    cout<<"Old = "<<losses<<"\n";
    
    // TODO: WE RESTRICT TO DELTA\in(-1,1) and Te in (1,20)
    double DeltaEff = Deltae;
    double TeEff  = Te;
    if(Deltae < -1.){
        DeltaEff = -1.;
    }
    if(Deltae > 1.){
        DeltaEff = 1.;
    }
    if(Te < 1.){
        TeEff = 1.;
    }
    if(TeEff > 20.){
        TeEff = 20.;
    }
    double Losses = 0.;
    double Delta_Delta = m_Delta[1] - m_Delta[0];
    double Delta_Temp  = m_Temperature[1] - m_Temperature[0];
    int lowerDeltaIdx = (int) (std::floor((DeltaEff - m_Delta[0])/Delta_Delta));
    int upperDeltaIdx = lowerDeltaIdx + 1;
    int lowerTempIdx = (int) (std::floor((TeEff-m_Temperature[0])/Delta_Temp));
    int upperTempIdx = lowerTempIdx + 1;

    double lowerDelta_lowerTemp = m_LossesFourthDelta[lowerDeltaIdx][ lowerTempIdx];
    double lowerDelta_upperTemp = m_LossesFourthDelta[lowerDeltaIdx][ upperTempIdx];
    double upperDelta_lowerTemp = m_LossesFourthDelta[upperDeltaIdx][ lowerTempIdx];
    double upperDelta_upperTemp = m_LossesFourthDelta[upperDeltaIdx][ upperTempIdx];
    
    double lower_Delta = linearInterpolation(m_Temperature[lowerTempIdx],m_Temperature[upperTempIdx] , lowerDelta_lowerTemp, lowerDelta_upperTemp, TeEff);
    double upper_Delta = linearInterpolation(m_Temperature[lowerTempIdx],m_Temperature[upperTempIdx] , upperDelta_lowerTemp, upperDelta_upperTemp, TeEff);
    Losses = linearInterpolation(m_Delta[lowerDeltaIdx],m_Delta[upperDeltaIdx] , lower_Delta, upper_Delta, DeltaEff);
    
    return -(1.)*Losses;
}

double SingleIonMixture::computeOmega11_en(const double& Te){

    double Omega11 = computeElecRate(Te, m_Omega11_en);
    
    return Omega11;
}

double SingleIonMixture::computeOmega12_en(const double& Te){

    double Omega12 = computeElecRate(Te, m_Omega12_en);
    
    return Omega12;
}

double SingleIonMixture::computeOmega13_en(const double& Te){

    double Omega13 = computeElecRate(Te, m_Omega13_en);
    
    return Omega13;
}

double SingleIonMixture::computeOmega14_en(const double& Te){

    double Omega14 = computeElecRate(Te, m_Omega14_en);
    
    return Omega14;
}

double SingleIonMixture::computeOmega22_ee(const double& Te){

    double Omega22 = computeElecRate(Te, m_Omega22_ee);
    
    return Omega22;
}

double SingleIonMixture::computeOmega23_ee(const double& Te){

    double Omega23 = computeElecRate(Te, m_Omega23_ee);
    
    return Omega23;
}

double SingleIonMixture::computeOmega24_ee(const double& Te){

    double Omega24 = computeElecRate(Te, m_Omega24_ee);
    
    return Omega24;
}

double SingleIonMixture::computeIonOmegaM(const double& T, const double& Mach){
    
    double Q_m = 0.;
    
    if(isnan(Mach)){
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        if(my_rank == MPI_WRITER){
            py::print("\n**************************************\n\n\n");
            py::print("\nMach number is NaN");
            py::print("\n**************************************\n\n\n");
        }
    }

    for (int i_Coll = 0; i_Coll < m_NumberIonColls; ++i_Coll){

        if(Mach >= m_IonMach.back()){
            int my_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            double last_Mach = m_IonMach.back();
            Q_m += computeIonRate(getOmegaMCoeffs(i_Coll, m_NumberMachEntries - 1), T);
            if(my_rank == MPI_WRITER){
                py::print("\n**************************************\n\n\n");
                py::print("\nWARNING: Fit for the IonOmegaM might be wrong as Mach = ",Mach," > ",last_Mach);
                py::print("\n**************************************\n\n\n");
            }
        }
        else if(Mach == 0.){
            Q_m += computeIonRate(getOmegaMCoeffs(i_Coll, 0), T);
        }
        else{
            int lowerIdx = (int) (std::floor(Mach/m_DeltaMach));
            int upperIdx = lowerIdx + 1;
            double lower_bound_Mach = getMach(i_Coll, lowerIdx);
            double upper_bound_Mach = getMach(i_Coll, upperIdx);
            double lowerBound = computeIonRate(getOmegaMCoeffs(i_Coll, lowerIdx), T);
            double upperBound = computeIonRate(getOmegaMCoeffs(i_Coll, upperIdx), T);

            Q_m += linearInterpolation(lower_bound_Mach, upper_bound_Mach, lowerBound, upperBound, Mach);
        }

    }

    return Q_m;
    
}

double SingleIonMixture::computeIonOmegaE(const double& T, const double& Mach){
    
    double Q_E = 0.;

    for (int i_Coll = 0; i_Coll < m_NumberIonColls; ++i_Coll){

        if(Mach >= m_IonMach.back()){
            int my_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            double last_Mach = m_IonMach.back();
            Q_E += computeIonRate(getOmegaECoeffs(i_Coll, m_NumberMachEntries - 1), T);
            if(my_rank == MPI_WRITER){
                py::print("\n**************************************\n\n\n");
                py::print("\nWARNING: Fit for the IonOmegaM might be wrong as Mach = ",Mach," > ",last_Mach);
                py::print("\n**************************************\n\n\n");
            }
        }
        else if(Mach == 0.){
            Q_E += computeIonRate(getOmegaECoeffs(i_Coll, 0), T);
        }
        else{
            int lowerIdx = (int) (std::floor(Mach/m_DeltaMach));
            int upperIdx = lowerIdx + 1;
            double lower_bound_Mach = getMach(i_Coll, lowerIdx);
            double upper_bound_Mach = getMach(i_Coll, upperIdx);
            double lowerBound = computeIonRate(getOmegaECoeffs(i_Coll, lowerIdx), T);
            double upperBound = computeIonRate(getOmegaECoeffs(i_Coll, upperIdx), T);

            Q_E += linearInterpolation(lower_bound_Mach, upper_bound_Mach, lowerBound, upperBound, Mach);
        }

    }
    return Q_E;
    
}
