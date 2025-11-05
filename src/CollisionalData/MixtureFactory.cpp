//
//  MixtureFactory.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#include "MixtureFactory.hpp"
#include <iostream>


using namespace std;

std::unique_ptr<Mixture> MixtureFactory::CreateMixture(const string &name)
{
    std::unique_ptr<Mixture> Mixture;
    if(name == "Null"){Mixture.reset(new NullMixture("NullMixture"));}
    if(name == "SingleIon"){Mixture.reset(new SingleIonMixture("SingleIon"));}
    else {
        py::print("Mixture Model not implemented. Taking No Mixture by default\n");
        Mixture.reset(new NullMixture("Null"));
    }
    
    return Mixture;
}

