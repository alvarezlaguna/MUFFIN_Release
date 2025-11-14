//
//  DataWriterFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "DataWriterFactory.hpp"

std::unique_ptr<DataWriter> DataWriterFactory::CreateDataWriter(const string &name)
{
    std::unique_ptr<DataWriter> dataWriter;
    if(name == "DataWriter1D"){dataWriter.reset(new DataWriter1D("DataWriter1D"));}
    else if(name == "DataWriter1DExtended"){dataWriter.reset(new DataWriter1DExtended("DataWriter1DExtended"));}
    else if(name == "DataWriter1DExtendedMPI"){dataWriter.reset(new DataWriter1DExtendedMPI("DataWriter1DExtendedMPI"));}
    else if(name == "DataWriter1DMPI"){dataWriter.reset(new DataWriter1DMPI("DataWriter1DMPI"));}
    else if(name == "DataWriter1DH5Py"){dataWriter.reset(new DataWriter1DH5Py("DataWriter1DH5Py"));}

    else {
        cout<<"Data Writer not implemented. Taking DataWriter1D by default\n";
        dataWriter.reset(new DataWriter1D("DataWriter1D"));
    }
    
    return dataWriter;
}

