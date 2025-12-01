//
//  FractionalTime2OSplitting.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FractionalTime2OSplitting_hpp
#define FractionalTime2OSplitting_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class FractionalTime2OSplitting : public TimeMethod
{
public:
    FractionalTime2OSplitting(string name) : TimeMethod(name), u_0(NBEQS*NBCELLS), u_1(NBEQS*NBCELLS) {}
    ~FractionalTime2OSplitting() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_0;
    vector<double> u_1;
};

#endif /* FractionalTime2OSplitting_hpp */
