//
//  NoReconstructor.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NoReconstructor_hpp
#define NoReconstructor_hpp

#include <stdio.h>
#include "SpaceReconstructor.hpp"
#include "../MeshData/MeshData.hpp"

class NoReconstructor : public SpaceReconstructor
{
public:
    NoReconstructor(string name) : SpaceReconstructor(name){}
    ~NoReconstructor() {}
    void setup(){}
    void unsetp(){}
    void reconstructField();
};

#endif /* NoReconstructor_hpp */
