//
//  OspreLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "OspreLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("ospre", OspreLimiter);

double OspreLimiter::operator()(const double& theta) {
    
    return 1.5*(theta*theta + theta)/(theta*theta + theta + 1);
}
