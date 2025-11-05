//
//  NullSourceTerm.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef NullSourceTerm_hpp
#define NullSourceTerm_hpp

#include <stdio.h>
#include <vector>
#include "SourceTerm.hpp"

using namespace Parameters;
using namespace std;

class NullSourceTerm : public SourceTerm
{
public:
    NullSourceTerm(string name) : SourceTerm(name) {}
    ~NullSourceTerm(){}
    virtual void computeSource();
};

#endif /* NullSourceTerm_hpp */
