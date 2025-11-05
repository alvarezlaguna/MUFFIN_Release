//
//  SSPRK5.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SSPRK5_hpp
#define SSPRK5_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class SSPRK5 : public TimeMethod
{
public:
    SSPRK5(string name) : TimeMethod(name), u_n(NBEQS*NBCELLS), u_1(NBEQS*NBCELLS), u_3(NBEQS*NBCELLS)  {}
    ~SSPRK5() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_n;
    vector<double> u_1;
    vector<double> u_3;
};

#endif /* SSPRK5_hpp */
