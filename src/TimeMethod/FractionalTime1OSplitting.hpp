//
//  FractionalTime1OSplitting.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FractionalTime1OSplitting_hpp
#define FractionalTime1OSplitting_hpp

#include <stdio.h>
#include <vector>
#include "TimeMethod.hpp"
#include "../Parameters.hpp"

using namespace std;
using namespace Parameters;


class FractionalTime1OSplitting : public TimeMethod
{
public:
    FractionalTime1OSplitting(string name) : TimeMethod(name), u_1(NBEQS*NBCELLS) {}
    ~FractionalTime1OSplitting() {}
    virtual void takeStep(double dt);
protected:
    vector<double> u_1;
};

#endif /* FractionalTime1OSplitting_hpp */
