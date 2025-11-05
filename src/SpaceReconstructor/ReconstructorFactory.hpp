//
//  ReconstructorFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ReconstructorFactory_hpp
#define ReconstructorFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>
#include "SpaceReconstructor.hpp"

using namespace std;

// Factory for creating instances of Reconstructor
class ReconstructorFactory
{
private:
    ReconstructorFactory();
    ~ReconstructorFactory() {}
    ReconstructorFactory(const ReconstructorFactory &) { }
    ReconstructorFactory &operator=(const ReconstructorFactory &) { return *this; }
    
    typedef std::function<std::unique_ptr<SpaceReconstructor>(const std::string&)> CreateReconstructorFn;
    std::map<std::string, CreateReconstructorFn> m_registry;
    std::mutex m_registryMutex;
public:
    //~ReconstructorFactory() { m_FactoryMap.clear(); }
    
    static ReconstructorFactory *Get()
    {
        static ReconstructorFactory instance;
        return &instance;
    }
    
    // Register a reconstructor creator function under the given name.
    // Returns true if registration succeeded (name was not already present).
    static bool Register(const std::string &name, CreateReconstructorFn pfnCreate);

    // Create a reconstructor by name. If name not found, returns a default NoReconstructor instance.
    static std::unique_ptr<SpaceReconstructor> CreateReconstructor(const std::string &name);
};

#endif /* ReconstructorFactory_hpp */
