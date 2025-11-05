//
//  PhysicalModelFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 05/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "PhysicalModelFactory.hpp"
#include "PhysicalModelRegistrar.hpp"
#include "EulerEq1D.hpp"
#include "PhysicalModelRegistrar.hpp"
using namespace std;

PhysicalModelFactory::PhysicalModelFactory() {}

bool PhysicalModelFactory::Register(const std::string &name, CreatePhysicalModelFn pfnCreate) {
    PhysicalModelFactory *instance = Get();
    std::lock_guard<std::mutex> lock(instance->m_registryMutex);
    bool success = instance->m_registry.emplace(name, pfnCreate).second;
    return success;
}

std::unique_ptr<PhysicalModel> PhysicalModelFactory::CreatePhysicalModel(const string &name)
{
    PhysicalModelFactory *instance = Get();
    std::lock_guard<std::mutex> lock(instance->m_registryMutex);
    
    auto it = instance->m_registry.find(name);
    if (it != instance->m_registry.end()) {
        return it->second(name);
    }
    
    cout << "Physical Model not registered: " << name << ". Taking Euler Eq. by default\n";
    return std::make_unique<EulerEq1D>("EulerEq1D");
}

