//
//  DiffusiveFluxFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef DiffusiveFluxFactory_hpp
#define DiffusiveFluxFactory_hpp

#include <stdio.h>
#include "DiffusiveFlux.hpp"
#include "NullDiffusiveFlux.hpp"
#include "FourierHeatFlux.hpp"


using namespace std;

// Factory for creating instances of Physical Models
class DiffusiveFluxFactory
{
private:
    DiffusiveFluxFactory();
    ~DiffusiveFluxFactory() {}
    DiffusiveFluxFactory(const DiffusiveFluxFactory &) { }
    DiffusiveFluxFactory &operator=(const DiffusiveFluxFactory &) { return *this; }
    
public:
    
    static DiffusiveFluxFactory *Get()
    {
        static DiffusiveFluxFactory instance;
        return &instance;
    }
    

    static std::unique_ptr<DiffusiveFlux> CreateDiffusiveFlux(const string &name);
};

#endif /* DiffusiveFluxFactory_hpp */
