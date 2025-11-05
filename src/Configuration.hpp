//
//  Configuration.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/03/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Configuration_hpp
#define Configuration_hpp

#include <stdio.h>
#include <pybind11/pybind11.h>
#include "Parameters.hpp"

namespace py = pybind11;

class Configuration{
    
public:
    Configuration(const py::dict &options);
};


#endif /* Configuration_hpp */
