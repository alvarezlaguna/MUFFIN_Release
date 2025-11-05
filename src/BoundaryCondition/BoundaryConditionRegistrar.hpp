#ifndef BoundaryConditionRegistrar_hpp
#define BoundaryConditionRegistrar_hpp

#include <iostream>
#include <functional>
#include "BoundaryCondition.hpp"
#include "BoundaryConditionFactory.hpp"

// Registration helper template
template<typename BCType>
class BoundaryConditionRegistrar {
public:
    BoundaryConditionRegistrar(const std::string& name) {
        auto creator = [](const std::string& bcName, const std::string& bcSide) -> std::unique_ptr<BoundaryCondition> {
            return std::make_unique<BCType>(bcName, bcSide);
        };
        
        bool success = BoundaryConditionFactory::Register(name, creator);
        if (success) {
            // std::cout << "Registered boundary condition: " << name << std::endl;
        } else {
            std::cerr << "Failed to register boundary condition: " << name << std::endl;
        }
    }
};

// Registration macro
#define REGISTER_BOUNDARY_CONDITION(name, classname) \
    static BoundaryConditionRegistrar<classname> g_register_##classname(name)

#endif // BoundaryConditionRegistrar_hpp