//
//  TVDReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVDReconstructor_hpp
#define TVDReconstructor_hpp

#include <stdio.h>
#include "SpaceReconstructor.hpp"
#include "Limiter/Limiter.hpp"
#include "Limiter/LimiterFactory.hpp"
#include <mpi.h>

class TVDReconstructor : public SpaceReconstructor
{
public:
    TVDReconstructor(string name) : SpaceReconstructor(name){}
    ~TVDReconstructor() {}
    void setup();
    void unsetup();
    virtual void reconstructField() = 0;
    
protected:
    std::unique_ptr<Limiter> m_limiter;
    
};

#endif /* TVDReconstructor_hpp */
