//
//  DataWriter1DExtended.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriter1DExtended_hpp
#define DataWriter1DExtended_hpp

#include <stdio.h>
#include "DataWriter.hpp"

class DataWriter1DExtended : public DataWriter
{
    
public:
    DataWriter1DExtended(string name) : DataWriter(name) {}
    ~DataWriter1DExtended() {}
    virtual void writeData(const int iter);
protected:
    unsigned int m_iter;
    
};


#endif /* DataWriter1DExtended_hpp */
