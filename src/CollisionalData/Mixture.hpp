//
//  Mixture.hpp
//  Fachade of a Mixture
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Mixture_hpp
#define Mixture_hpp

#include <stdio.h>
#include <vector>
#include "../Component.hpp"

using namespace Parameters;
using namespace std;

class Mixture : public Component
{
public :
    Mixture(string name) : Component(name) {}
    virtual ~Mixture(){}
    
    virtual void setup() {};
    virtual void unsetup() {};
};
#endif /* Mixture_hpp */
