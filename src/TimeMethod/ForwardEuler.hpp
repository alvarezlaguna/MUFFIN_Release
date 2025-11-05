//
//  ForwardEuler.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ForwardEuler_hpp
#define ForwardEuler_hpp

#include <stdio.h>
#include "TimeMethod.hpp"


class ForwardEuler : public TimeMethod
{
public:
    ForwardEuler(string name) : TimeMethod(name){}
    ~ForwardEuler() {}
    virtual void takeStep(double dt);

};

#endif /* ForwardEuler_hpp */
