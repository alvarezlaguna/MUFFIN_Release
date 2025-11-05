//
//  NullMixture.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NullMixture_hpp
#define NullMixture_hpp

#include <stdio.h>
#include "Mixture.hpp"

using namespace Parameters;

class NullMixture : public Mixture
{
public:
    NullMixture(string name) : Mixture(name) {}
    ~NullMixture(){}
    
protected:

};

#endif /* NullMixture_hpp */
