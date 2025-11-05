//
//  TVD2ndOrder1DReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVD2ndOrder1DReconstructor_hpp
#define TVD2ndOrder1DReconstructor_hpp

#include <stdio.h>
#include "TVDReconstructor.hpp"
#include "../MeshData/MeshData.hpp"

class TVD2ndOrder1DReconstructor : public TVDReconstructor
{
public:
    TVD2ndOrder1DReconstructor(string name): TVDReconstructor(name){}
    ~TVD2ndOrder1DReconstructor() {}
    void setup(){TVDReconstructor::setup();};
    void unsetup(){TVDReconstructor::unsetup();}
    virtual void reconstructField();
};

#endif /* TVD2ndOrder1DReconstructor_hpp */
