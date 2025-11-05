//
//  SSPRK3.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SSPRK3_hpp
#define SSPRK3_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class SSPRK3 : public TimeMethod
{
public:
    SSPRK3(string name) : TimeMethod(name), u_1(NBEQS*NBCELLS) {}
    ~SSPRK3() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_1;
};

#endif /* SSPRK3_hpp */
