//
//  BoundaryConditionFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef BoundaryConditionFactory_hpp
#define BoundaryConditionFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <mutex>
#include <memory>
#include <functional>
#include "BoundaryCondition.hpp"

using namespace std;

// Factory for creating instances of Boundary Conditions
class BoundaryConditionFactory {
public:
    using Creator = std::function<std::unique_ptr<BoundaryCondition>(const std::string&, const std::string&)>;

    static BoundaryConditionFactory& getInstance();

    // Static registration method
    static bool Register(const std::string& name, Creator creator);

    // Create method
    static std::unique_ptr<BoundaryCondition> CreateBoundaryCondition(const std::string& name, const std::string& side);

private:
    BoundaryConditionFactory() = default;

    // Registry of boundary condition creators
    std::map<std::string, Creator> m_registry;
    std::mutex m_mutex;

    // Singleton instance
    static BoundaryConditionFactory* m_instance;
    static std::mutex m_instanceMutex;
    BoundaryConditionFactory(const BoundaryConditionFactory &) { }
    BoundaryConditionFactory &operator=(const BoundaryConditionFactory &) { return *this; }

public:
    static BoundaryConditionFactory *Get()
    {
        static BoundaryConditionFactory instance;
        return &instance;
    }
};

#endif /* BoundaryConditionFactory_hpp */
