//
//  SpaceReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SpaceReconstructor_hpp
#define SpaceReconstructor_hpp

#include <stdio.h>
#include "../Component.hpp"

class SpaceReconstructor : public Component
{
public:
    SpaceReconstructor(string name) : Component(name) {}
    virtual ~SpaceReconstructor(){}
    void setup(){}
    void unsetup(){}
    virtual void reconstructField() = 0;
};

#endif /* SpaceReconstructor_hpp */
