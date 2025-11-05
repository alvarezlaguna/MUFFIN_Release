//
//  vanAlbadaLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "vanAlbadaLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("vanAlbada", vanAlbadaLimiter);
#include <cmath>

using namespace std;

double vanAlbadaLimiter::operator()(const double& theta) {
    
    return (theta*theta + theta)/(theta*theta + 1); // van Aldaba 2 limiter
    //return 1.5*(theta + abs(theta))/(theta + 2); //HCUS
    //return (theta + abs(theta))/(abs(theta) + 1); // van Aldaba 2 limiter
}
