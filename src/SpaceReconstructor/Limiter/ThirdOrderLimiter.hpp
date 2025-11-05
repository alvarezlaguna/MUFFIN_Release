//
//  ThirdOrderLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ThirdOrderLimiter_hpp
#define ThirdOrderLimiter_hpp

#include <stdio.h>
#include "Limiter.hpp"

class ThirdOrderLimiter : public Limiter
{
public:
    ThirdOrderLimiter(string name) : Limiter(name){}
    ~ThirdOrderLimiter() {}
    virtual double operator()(const double& theta);
};

#endif /* ThirdOrderLimiter_hpp */
