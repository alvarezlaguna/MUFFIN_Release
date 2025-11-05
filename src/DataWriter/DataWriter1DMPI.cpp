//
//  DataWriter1DMPI.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/04/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DataWriter1DMPI.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../MeshData/MeshData.hpp"
#include "../Parameters.hpp"
#include "../MeshData/Cell1D.hpp"

#include <mpi.h>


using namespace Parameters;
using namespace std;

void DataWriter1DMPI::writeData(const int iter)
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    vector<Cell1D>& cells  = MeshData::getInstance().getData<Cell1D>("Cells");
    vector<double>& x = MeshData::getInstance().getData<double>("x");
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    // Write the data in ASCII
    ofstream dataFile;      // file object
    char fname[200];        //buffer for the file name
    m_iter = iter;
    
    int rank, i; char line[200];
    
    MPI_File fh ;
    MPI_Status status ;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Write the file name
    sprintf(fname, "%s", RESULTDIRECTORY.c_str());
    sprintf(fname + strlen(fname), "/file_iter_%06i_time_%.4e.txt", m_iter, physTime); //TODO writing the name in a general way
    
    // Open the file with MPI
    MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    
    int err, count;
    
    for (int iCell = 0; iCell < NBCELLS; ++iCell){
        int offset = 0;
        for (int iP = 0; iP < rank; ++iP){
            offset += NBCELLS_MPI[iP];
        }
        offset += iCell;
        
        count = sprintf(line, "%*.8e\t",18, x[iCell]);
        
        int sizeOfLine = count*(NBEQS + 1) + 1; // The full line is the x plus the variables plus the newline character
        // Write the x
        MPI_Offset displace = offset*sizeOfLine*sizeof(char); // start of the view for each processor
        err = MPI_File_write_at(fh,displace, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
        for(int iEq = 0; iEq < NBEQS ; ++iEq ){
            // Write the variables
            count = sprintf(line, "%*.8e\t", 18, cells[iCell].uCC[iEq]); // Line printing the result
            displace += count*sizeof(char); // start of the view for each processor
            err = MPI_File_write_at(fh,displace, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        
        // Write en line
        displace += count*sizeof(char); // start of the view for each processor
        count = sprintf(line, "\n");                             // End line
        err = MPI_File_write_at(fh,displace, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    MPI_File_close(&fh ) ;
    
}

