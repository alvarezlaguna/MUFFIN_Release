//
//  MixtureFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef MixtureFactory_hpp
#define MixtureFactory_hpp

#include <stdio.h>
#include "Mixture.hpp"
#include "SingleIonMixture.hpp"
#include "NullMixture.hpp"


using namespace std;

// Factory for creating instances of Physical Models
class MixtureFactory
{
private:
    MixtureFactory();
    ~MixtureFactory() {}
    MixtureFactory(const MixtureFactory &) { }
    MixtureFactory &operator=(const MixtureFactory &) { return *this; }
    
public:
    
    static MixtureFactory *Get()
    {
        static MixtureFactory instance;
        return &instance;
    }
    

    static std::unique_ptr<Mixture> CreateMixture(const string &name);
};

#endif /* MixtureFactory_hpp */
