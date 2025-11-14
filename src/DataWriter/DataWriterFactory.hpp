//
//  DataWriterFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DataWriterFactory_hpp
#define DataWriterFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include "DataWriter.hpp"
#include "DataWriter1D.hpp"
#include "DataWriter1DExtended.hpp"
#include "DataWriter1DExtendedMPI.hpp"
#include "DataWriter1DMPI.hpp"
#include "DataWriter1DH5Py.hpp"

using namespace std;

// Factory for creating instances of Limiter
class DataWriterFactory
{
private:
    DataWriterFactory();
    DataWriterFactory(const DataWriterFactory &) { }
    DataWriterFactory &operator=(const DataWriterFactory &) { return *this; }
    
    //typedef map<string, void*> FactoryMap;
    //FactoryMap m_FactoryMap;
public:
    //~LimiterFactory() { m_FactoryMap.clear(); }
    
    static DataWriterFactory *Get()
    {
        static DataWriterFactory instance;
        return &instance;
    }
    
    //void Register(const string &name, CreateLimiterFn pfnCreate);
    static std::unique_ptr<DataWriter> CreateDataWriter(const string &name);
};

#endif /* DataWriterFactory_hpp */
