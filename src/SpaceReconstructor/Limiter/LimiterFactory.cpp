//
//  LimiterFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "LimiterFactory.hpp"
using namespace std;
// Provide definition of private constructor
LimiterFactory::LimiterFactory() {}

#include "OspreLimiter.hpp" // default fallback

bool LimiterFactory::Register(const std::string &name, CreateLimiterFn pfnCreate)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    if(mapRef.find(name) != mapRef.end()) {
        std::cout << "LimiterFactory: registration FAILED - '" << name << "' is already registered\n";
        return false;
    }
    mapRef[name] = pfnCreate;
    // std::cout << "LimiterFactory: registered limiter '" << name << "' successfully\n";
    return true;
}

std::unique_ptr<Limiter> LimiterFactory::CreateLimiter(const string &name)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    auto it = mapRef.find(name);
    if(it != mapRef.end()){
        return (it->second)(name);
    }

    std::cout << "Limiter not implemented: '" << name << "'. Taking 'ospre' by default\n";
    return std::unique_ptr<Limiter>(new OspreLimiter("ospre"));
}