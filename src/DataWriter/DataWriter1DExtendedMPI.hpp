//
//  DataWriter1DExtendedMPI.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriter1DExtendedMPI_hpp
#define DataWriter1DExtendedMPI_hpp

#include <stdio.h>
#include "DataWriter.hpp"

class DataWriter1DExtendedMPI : public DataWriter
{
    
public:
    DataWriter1DExtendedMPI(string name) : DataWriter(name) {}
    ~DataWriter1DExtendedMPI() {}
    virtual void writeData(const int iter);
protected:
    unsigned int m_iter;
    
};
#endif /* DataWriter1DExtendedMPI_hpp */
