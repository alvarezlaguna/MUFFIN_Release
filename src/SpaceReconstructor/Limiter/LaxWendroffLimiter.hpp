//
//  LaxWendroffLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LaxWendroffLimiter_hpp
#define LaxWendroffLimiter_hpp

#include <stdio.h>
#include "Limiter.hpp"

class LaxWendroffLimiter : public Limiter
{
public:
    LaxWendroffLimiter(string name) : Limiter(name){}
    ~LaxWendroffLimiter() {}
    virtual double operator()(const double& theta);
};


#endif /* LaxWendroffLimiter_hpp */
