//
//  DataWriter1DExtended.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DataWriter1DExtended.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../MeshData/MeshData.hpp"
#include "../Parameters.hpp"
#include "../MeshData/Cell1D.hpp"


using namespace Parameters;
using namespace std;

void DataWriter1DExtended::writeData(const int iter)
{
    vector<Cell1D>& cells  = MeshData::getInstance().getData<Cell1D>("Cells");
    vector<double>& phi  = MeshData::getInstance().getData<double>("Phi");
    vector<double>& x = MeshData::getInstance().getData<double>("x");
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
//    vector<double>& nu_bsc_ionXe  = MeshData::getInstance().getData<double>("nu_bsc_ionXe");
//    vector<double>& nu_ela_elecXe = MeshData::getInstance().getData<double>("nu_ela_elecXe");
//    vector<double>& nu_exc_elecXe = MeshData::getInstance().getData<double>("nu_exc_elecXe");
//    vector<double>& nu_iso_ionXe  = MeshData::getInstance().getData<double>("nu_iso_ionXe");
//    vector<double>& nu_iz_elecXe  = MeshData::getInstance().getData<double>("nu_iz_elecXe");
//
//    vector<double>& R_xx          = MeshData::getInstance().getData<double>("R_xx");
//
//    vector<double>& Jez           = MeshData::getInstance().getData<double>("Jez");
//
    // Write the data in ASCII
    ofstream dataFile;      // file object
    char fname[200];        //buffer
    m_iter = iter;
    sprintf(fname, "%s", RESULTDIRECTORY.c_str());
    sprintf(fname + strlen(fname), "/file_iter_%06i_time_%.8e.txt", m_iter, physTime); //TODO writing the name in a general way
    dataFile.open(fname);
    for (unsigned int i = 0; i < NBCELLS; i++)
    {
        dataFile << x[i] << "\t";
        for(unsigned int iEq = 0; iEq < NBEQS; iEq++){
            dataFile << cells[i].u_CC()[iEq]<<"\t";
        }
        dataFile << phi[i]<<"\t";
//        dataFile << phi[i]<<"\t";
//        dataFile << nu_bsc_ionXe[i]<<"\t";
//        dataFile << nu_iso_ionXe[i]<<"\t";
//        dataFile << nu_ela_elecXe[i]<<"\t";
//        dataFile << nu_exc_elecXe[i]<<"\t";
//        dataFile << nu_iz_elecXe[i]<<"\t";
//        dataFile << R_xx[i]<<"\t";
//        dataFile << Jez[i]<<"\t";
        dataFile <<"\n";
        
    }
    dataFile.close();
}
