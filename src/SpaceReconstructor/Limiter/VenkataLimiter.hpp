//
//  VenkataLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef VenkataLimiter_hpp
#define VenkataLimiter_hpp

#include <stdio.h>
#include "Limiter.hpp"

class VenkataLimiter : public Limiter
{
public:
    VenkataLimiter(string name) : Limiter(name){}
    ~VenkataLimiter() {}
    virtual double operator()(const double& theta);
};

#endif /* VenkataLimiter_hpp */
