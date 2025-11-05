//
//  SourceTermFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "SourceTermFactory.hpp"
// Define private constructor
SourceTermFactory::SourceTermFactory() {}

#include "NullSourceTerm.hpp" // default fallback

bool SourceTermFactory::Register(const std::string &name, CreateSourceTermFn pfnCreate)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    if(mapRef.find(name) != mapRef.end()){
        std::cout << "SourceTermFactory: registration FAILED - '" << name << "' is already registered\n";
        return false;
    }
    mapRef[name] = pfnCreate;
    // std::cout << "SourceTermFactory: registered source-term '" << name << "' successfully\n";
    return true;
}

std::unique_ptr<SourceTerm> SourceTermFactory::CreateSourceTerm(const string &name)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    auto it = mapRef.find(name);
    if(it != mapRef.end()){
        return (it->second)(name);
    }

    std::cout << "Source Term not implemented: '" << name << "'. Taking 'NullSourceTerm' by default\n";
    return std::unique_ptr<SourceTerm>(new NullSourceTerm("NullSourceTerm"));
}
