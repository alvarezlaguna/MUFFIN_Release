//
//  StrangSplitting.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef StrangSplitting_hpp
#define StrangSplitting_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class StrangSplitting : public TimeMethod
{
public:
    StrangSplitting(string name) : TimeMethod(name), u_0(NBEQS*NBCELLS), u_1(NBEQS*NBCELLS) {}
    ~StrangSplitting() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_0;
    vector<double> u_1;
};

#endif /* StrangSplitting_hpp */
