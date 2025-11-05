//
//  OspreLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef OspreLimiter_hpp
#define OspreLimiter_hpp

#include <stdio.h>
#include "Limiter.hpp"

class OspreLimiter : public Limiter
{
public:
    OspreLimiter(string name) : Limiter(name){}
    ~OspreLimiter() {}
    virtual double operator()(const double& theta);
};


#endif /* OspreLimiter_hpp */
