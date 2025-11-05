//
//  VenkataLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "VenkataLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("Venkatakrishnan", VenkataLimiter);

double VenkataLimiter::operator()(const double& theta) {
    
    return (theta*theta + 2*theta)/(theta*theta + theta + 2); // Venkatakrishnan limiter
}