//
//  ReconstructorFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ReconstructorFactory.hpp"
using namespace std;
// Define the (previously only-declared) private constructor so the
// singleton instance can be created. This resolves the linker error for
// symbol '__ZN20ReconstructorFactoryC1Ev' (ReconstructorFactory::ReconstructorFactory()).
ReconstructorFactory::ReconstructorFactory() {}

#include "NoReconstructor.hpp" // for default fallback

bool ReconstructorFactory::Register(const std::string &name, CreateReconstructorFn pfnCreate)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    if(mapRef.find(name) != mapRef.end()) {
        std::cout << "ReconstructorFactory: registration FAILED - '" << name << "' is already registered\n";
        return false;
    }
    mapRef[name] = pfnCreate;
    //std::cout << "ReconstructorFactory: registered reconstructor '" << name << "' successfully\n";
    return true;
}

std::unique_ptr<SpaceReconstructor> ReconstructorFactory::CreateReconstructor(const string &name)
{
    std::lock_guard<std::mutex> lock(Get()->m_registryMutex);
    auto &mapRef = Get()->m_registry;
    auto it = mapRef.find(name);
    if(it != mapRef.end()){
        return (it->second)(name);
    }

    std::cout << "Space reconstruction method '" << name << "' not implemented. Taking no reconstructor by default\n";
    return std::unique_ptr<SpaceReconstructor>(new NoReconstructor("1stOrder"));
}
