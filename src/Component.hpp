//
//  Component.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Component_hpp
#define Component_hpp

#include <stdio.h>
#include "Parameters.hpp"
#include <string>
#include <pybind11/pybind11.h>


namespace py = pybind11;
using namespace std;
using namespace Parameters;

class Component
{
public :
    Component(string name) : s(name) {} //constructor
    virtual ~Component() {}
    virtual void setup() {}
    virtual void unsetup() {}
    string getName() const { return s;}
    
protected:
    string s;

};

#endif /* Component_hpp */
