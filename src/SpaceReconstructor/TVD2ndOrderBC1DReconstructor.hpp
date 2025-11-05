//
//  TVD2ndOrderBC1DReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 04/12/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVD2ndOrderBC1DReconstructor_hpp
#define TVD2ndOrderBC1DReconstructor_hpp

#include <stdio.h>
#include "TVDReconstructor.hpp"
#include "../MeshData/MeshData.hpp"

class TVD2ndOrderBC1DReconstructor : public TVDReconstructor
{
public:
    TVD2ndOrderBC1DReconstructor(string name): TVDReconstructor(name){}
    ~TVD2ndOrderBC1DReconstructor() {}
    void setup(){TVDReconstructor::setup();};
    void unsetup(){TVDReconstructor::unsetup();}
    virtual void reconstructField();
};

#endif /* TVD2ndOrderBC1DReconstructor_hpp */
