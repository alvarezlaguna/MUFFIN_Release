//
//  NullDiffusiveFlux.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NullDiffusiveFlux_hpp
#define NullDiffusiveFlux_hpp

#include <stdio.h>
#include "DiffusiveFlux.hpp"

using namespace Parameters;

class NullDiffusiveFlux : public DiffusiveFlux
{
public:
    NullDiffusiveFlux(string name) : DiffusiveFlux(name) {}
    ~NullDiffusiveFlux(){}
    
protected:

};

#endif /* NullDiffusiveFlux_hpp */
