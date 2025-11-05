//
//  LaxWendroffLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 27/11/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "LaxWendroffLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("LaxWendroff", LaxWendroffLimiter);

double LaxWendroffLimiter::operator()(const double& theta) {
    
    return 1.;
}
