//
//  TVDRK3.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TVDRK3_hpp
#define TVDRK3_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class TVDRK3 : public TimeMethod
{
public:
    TVDRK3(string name) : TimeMethod(name), u_1(NBEQS*NBCELLS) {}
    ~TVDRK3() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_1;
};

#endif /* TVDRK3_hpp */
