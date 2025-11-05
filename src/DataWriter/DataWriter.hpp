//
//  DataWriter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriter_hpp
#define DataWriter_hpp

#include <stdio.h>
#include "../Component.hpp"

class DataWriter : public Component
{
public:
    DataWriter(string name) : Component(name) {}
    virtual ~DataWriter(){}
    virtual void writeData(const int iter) = 0;
    
protected:
    unsigned int m_iter;

};

#endif /* DataWriter_hpp */
