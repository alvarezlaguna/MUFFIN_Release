//
//  Limiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Limiter_hpp
#define Limiter_hpp

#include <stdio.h>
#include "../../Component.hpp"

class Limiter : public Component
{
public:
    Limiter(string name) : Component(name) {}
    virtual ~Limiter() {}
    virtual double operator()(const double& theta) = 0;
};

#endif /* Limiter_hpp */
