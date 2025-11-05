//
//  FluxSchemeRegistrar.hpp
//  Muffin
//
//  Created on 2025-11-03.
//

#ifndef FluxSchemeRegistrar_hpp
#define FluxSchemeRegistrar_hpp

#include <stdio.h>
#include <iostream>
#include <memory>
#include "FluxSchemeFactory.hpp"

// Helper class for flux scheme registration
template<typename T>
class FluxSchemeRegistrar {
public:
    FluxSchemeRegistrar(const std::string &name) {
        bool registered = FluxSchemeFactory::Register(name,
            [](const std::string &n) -> std::unique_ptr<FluxScheme> {
                return std::make_unique<T>(n);
            }
        );
        if (registered) {
            // std::cout << "Successfully registered flux scheme: " << name << std::endl;
        } else {
            std::cout << "Warning: Failed to register flux scheme: " << name << " (already exists)" << std::endl;
        }
    }
};

// Helper macro to register a flux scheme
#define REGISTER_FLUXSCHEME(name, classname) \
    static FluxSchemeRegistrar<classname> g_register_##classname(name)

#endif /* FluxSchemeRegistrar_hpp */
