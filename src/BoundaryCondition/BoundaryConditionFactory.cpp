//
//  BoundaryConditionFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "BoundaryConditionFactory.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include "NeumannBC.hpp"
#include "BoundaryConditionRegistrar.hpp"
#include <stdexcept>

// Initialize static members
BoundaryConditionFactory* BoundaryConditionFactory::m_instance = nullptr;
std::mutex BoundaryConditionFactory::m_instanceMutex;

BoundaryConditionFactory& BoundaryConditionFactory::getInstance() {
    std::lock_guard<std::mutex> lock(m_instanceMutex);
    if (m_instance == nullptr) {
        m_instance = new BoundaryConditionFactory();
    }
    return *m_instance;
}

bool BoundaryConditionFactory::Register(const std::string& name, Creator creator) {
    auto& instance = getInstance();
    std::lock_guard<std::mutex> lock(instance.m_mutex);
    
    if (instance.m_registry.find(name) != instance.m_registry.end()) {
        std::cerr << "Warning: Boundary Condition '" << name << "' already registered" << std::endl;
        return false;
    }
    
    instance.m_registry[name] = creator;
    return true;
}

std::unique_ptr<BoundaryCondition> BoundaryConditionFactory::CreateBoundaryCondition(const std::string& name, const std::string& side) {
    auto& instance = getInstance();
    std::lock_guard<std::mutex> lock(instance.m_mutex);
    
    auto it = instance.m_registry.find(name);
    if (it != instance.m_registry.end()) {
        return it->second(name, side);
    }
    
    std::cerr << "Warning: Unknown boundary condition '" << name 
              << "'. Using NeumannBC." << std::endl;
              
    // Return default Neumann boundary condition
    return std::make_unique<NeumannBC>(name, side);
}
