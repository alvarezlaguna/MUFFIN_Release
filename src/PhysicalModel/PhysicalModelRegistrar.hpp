//
//  PhysicalModelRegistrar.hpp
//  Muffin
//
//  Created by GitHub Copilot on 02/11/25.
//  Copyright © 2025 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PhysicalModelRegistrar_hpp
#define PhysicalModelRegistrar_hpp

#include <stdio.h>
#include <iostream>
#include <memory>
#include "PhysicalModelFactory.hpp"

// Helper class for physical model registration
template<typename T>
class PhysicalModelRegistrar {
public:
    PhysicalModelRegistrar(const std::string &name) {
        bool registered = PhysicalModelFactory::Register(name,
            [](const std::string &n) -> std::unique_ptr<PhysicalModel> {
                return std::make_unique<T>(n);
            }
        );
        if (registered) {
            // std::cout << "Successfully registered physical model: " << name << std::endl;
        } else {
            std::cout << "Warning: Failed to register physical model: " << name << " (already exists)" << std::endl;
        }
    }
};

// Helper macro to register a physical model
#define REGISTER_PHYSICALMODEL(name, classname) \
    static PhysicalModelRegistrar<classname> g_register_##classname(name)

#endif /* PhysicalModelRegistrar_hpp */