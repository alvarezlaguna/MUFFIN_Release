//
//  DataWriter1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriter1D_hpp
#define DataWriter1D_hpp

#include <stdio.h>
#include "DataWriter.hpp"

class DataWriter1D : public DataWriter
{
    
public:
    DataWriter1D(string name) : DataWriter(name) {}
    ~DataWriter1D() {}
    virtual void writeData(const int iter);
protected:
    unsigned int m_iter;
    
};


#endif /* DataWriter1D_hpp */
