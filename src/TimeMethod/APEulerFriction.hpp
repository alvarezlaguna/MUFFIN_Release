//
//  APEulerFriction.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef APEulerFriction_hpp
#define APEulerFriction_hpp

#include <stdio.h>
#include "TimeMethod.hpp"


class APEulerFriction : public TimeMethod
{
public:
    APEulerFriction(string name) : TimeMethod(name){}
    ~APEulerFriction() {}
    virtual void takeStep(double dt);
};


#endif /* APEulerFriction_hpp */
