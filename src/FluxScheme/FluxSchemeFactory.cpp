//
//  FluxSchemeFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FluxSchemeFactory.hpp"
#include "LaxFriedrichFlux.hpp"
#include <mutex>

FluxSchemeFactory::FluxSchemeFactory() {}

bool FluxSchemeFactory::Register(const std::string &name, CreateFluxSchemeFn pfnCreate) {
    FluxSchemeFactory *instance = Get();
    std::lock_guard<std::mutex> lock(instance->m_registryMutex);
    bool success = instance->m_registry.emplace(name, pfnCreate).second;
    return success;
}

std::unique_ptr<FluxScheme> FluxSchemeFactory::CreateFluxScheme(const string &name)
{
    FluxSchemeFactory *instance = Get();
    {
        std::lock_guard<std::mutex> lock(instance->m_registryMutex);
        auto it = instance->m_registry.find(name);
        if (it != instance->m_registry.end()) {
            return it->second(name);
        }
    }

    cout << "Flux Scheme not registered: " << name << ". Taking LaxFriedrich by default\n";
    return std::make_unique<LaxFriedrichFlux>("LaxFriedrich");

}
