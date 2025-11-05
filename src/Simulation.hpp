//
//  Simulation.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Simulation_hpp
#define Simulation_hpp

#include <stdio.h>
#include "Component.hpp"
#include "FluxScheme/FluxScheme.hpp"
#include "SpaceReconstructor/SpaceReconstructor.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"


class Simulation : public Component {
public:
    Simulation(string name) : Component(name) {}
    virtual ~Simulation() {}
    void setup();
    void unsetup();
    virtual pybind11::array_t<double> simulate() = 0;
};

#endif /* Simulation_hpp */
