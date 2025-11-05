//
//  FluxSchemeFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FluxSchemeFactory_hpp
#define FluxSchemeFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>
#include "FluxScheme.hpp"


using namespace std;

// Factory for creating instances of FluxScheme
class FluxSchemeFactory
{
private:
    FluxSchemeFactory();
    ~FluxSchemeFactory() {}
    FluxSchemeFactory(const FluxSchemeFactory &) { }
    FluxSchemeFactory &operator=(const FluxSchemeFactory &) { return *this; }
    
    typedef std::function<std::unique_ptr<FluxScheme>(const std::string&)> CreateFluxSchemeFn;
    std::map<std::string, CreateFluxSchemeFn> m_registry;
    std::mutex m_registryMutex;
public:
    //~LimiterFactory() { m_FactoryMap.clear(); }
    
    static FluxSchemeFactory *Get()
    {
        static FluxSchemeFactory instance;
        return &instance;
    }
    
    // Register a flux-scheme creator function under the given name.
    // Returns true if registration succeeded (name was not already present).
    static bool Register(const std::string &name, CreateFluxSchemeFn pfnCreate);

    // Create a flux scheme by name. If name not found, returns a default LaxFriedrich instance.
    static std::unique_ptr<FluxScheme> CreateFluxScheme(const std::string &name);
};

#endif /* FluxSchemeFactory_hpp */
