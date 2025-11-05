//
//  vanAlbadaLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef vanAlbadaLimiter_hpp
#define vanAlbadaLimiter_hpp

#include <stdio.h>
#include "Limiter.hpp"

class vanAlbadaLimiter : public Limiter
{
public:
    vanAlbadaLimiter(string name) : Limiter(name){}
    ~vanAlbadaLimiter() {}
    virtual double operator()(const double& theta);
};

#endif /* vanAlbadaLimiter_hpp */
