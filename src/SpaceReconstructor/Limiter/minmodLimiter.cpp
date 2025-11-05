//
//  minmodLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "minmodLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("minmod", minmodLimiter);
#include <cmath>
#include <algorithm>   

double minmodLimiter::operator()(const double& theta) {
    
    return std::max(0., std::min(1.,theta)); // minmod limiter
}
