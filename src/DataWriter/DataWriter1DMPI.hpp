//
//  DataWriter1DMPI.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/04/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriter1DMPI_hpp
#define DataWriter1DMPI_hpp

#include <stdio.h>
#include "DataWriter.hpp"

class DataWriter1DMPI : public DataWriter
{
    
public:
    DataWriter1DMPI(string name) : DataWriter(name) {}
    ~DataWriter1DMPI() {}
    virtual void writeData(const int iter);
protected:
    unsigned int m_iter;
    
};
#endif /* DataWriter1DExtendedMPI_hpp */
