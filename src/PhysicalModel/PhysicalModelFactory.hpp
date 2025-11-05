//
//  PhysicalModelFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 05/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PhysicalModelFactory_hpp
#define PhysicalModelFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>
#include "PhysicalModel.hpp"

using namespace std;

// Factory for creating instances of Physical Models
class PhysicalModelFactory
{
private:
    PhysicalModelFactory();
    ~PhysicalModelFactory() {}
    PhysicalModelFactory(const PhysicalModelFactory &) { }
    PhysicalModelFactory &operator=(const PhysicalModelFactory &) { return *this; }
    
    typedef std::function<std::unique_ptr<PhysicalModel>(const std::string&)> CreatePhysicalModelFn;
    std::map<std::string, CreatePhysicalModelFn> m_registry;
    std::mutex m_registryMutex;
    
public:
    static PhysicalModelFactory *Get()
    {
        static PhysicalModelFactory instance;
        return &instance;
    }
    
    // Register a physical model creator function under the given name.
    // Returns true if registration succeeded (name was not already present).
    static bool Register(const std::string &name, CreatePhysicalModelFn pfnCreate);

    // Create a physical model by name. If name not found, returns an EulerEq1D instance.
    static std::unique_ptr<PhysicalModel> CreatePhysicalModel(const std::string &name);
};

#endif /* PhysicalModelFactory_hpp */
