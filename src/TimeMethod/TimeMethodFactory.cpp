//
//  TimeMethodFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TimeMethodFactory.hpp"
#include "ForwardEuler.hpp" // default fallback

// Define private constructor
TimeMethodFactory::TimeMethodFactory() {}

bool TimeMethodFactory::Register(const std::string &name, CreateTimeMethodFn pfnCreate)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    if(mapRef.find(name) != mapRef.end()){
        std::cout << "TimeMethodFactory: registration FAILED - '" << name << "' is already registered\n";
        return false;
    }
    mapRef[name] = pfnCreate;
    // std::cout << "TimeMethodFactory: registered time method '" << name << "' successfully\n";
    return true;
}

std::unique_ptr<TimeMethod> TimeMethodFactory::CreateTimeMethod(const string &name)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    auto it = mapRef.find(name);
    if(it != mapRef.end()){
        return (it->second)(name);
    }

    std::cout << "Time Method not implemented: '" << name << "'. Taking 'ForwardEuler' by default\n";
    return std::unique_ptr<TimeMethod>(new ForwardEuler("ForwardEuler"));
}
