//
//  DataWriter1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>

#include "DataWriter1D.hpp"
#include "../MeshData/MeshData.hpp"
#include "../Parameters.hpp"
#include "../MeshData/Cell1D.hpp"


using namespace Parameters;
using namespace std;

void DataWriter1D::writeData(const int iter)
{
    vector<Cell1D>& cells  = MeshData::getInstance().getData<Cell1D>("Cells");
    //vector<double>& phi  = MeshData::getInstance().getData<double>("Phi");
    vector<double>& x = MeshData::getInstance().getData<double>("x");
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    // Write the data in ASCII
    ofstream dataFile;      // file object
    char fname[200];        //buffer
    m_iter = iter;
    sprintf(fname, "%s", RESULTDIRECTORY.c_str());
    sprintf(fname + strlen(fname), "/file_iter_%06i_time_%.4e.txt", m_iter, physTime); //TODO writing the name in a general way
    dataFile.open(fname);
    for (unsigned int i = 0; i < NBCELLS; i++)
    {
        dataFile << x[i] << "\t";
        for(unsigned int iEq = 0; iEq < NBEQS; iEq++){
            dataFile << cells[i].u_CC()[iEq]<<"\t";
        }
        //dataFile << phi[i]<<"\t";
        dataFile <<"\n";
        
    }
    dataFile.close();
}
