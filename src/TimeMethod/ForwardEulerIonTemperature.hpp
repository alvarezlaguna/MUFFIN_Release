//
//  ForwardEulerIonTemperature.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ForwardEulerIonTemperature_hpp
#define ForwardEulerIonTemperature_hpp

#include <stdio.h>
#include "TimeMethod.hpp"


class ForwardEulerIonTemperature : public TimeMethod
{
public:
    ForwardEulerIonTemperature(string name) : TimeMethod(name){}
    ~ForwardEulerIonTemperature() {}
    virtual void takeStep(double dt);
};

#endif /* ForwardEulerIonTemperature_hpp */
