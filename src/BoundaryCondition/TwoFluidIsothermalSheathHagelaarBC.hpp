//
//  TwoFluidIsothermalSheathHagelaarBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermalSheathHagelaarBC_hpp
#define TwoFluidIsothermalSheathHagelaarBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class TwoFluidIsothermalSheathHagelaarBC : public BoundaryCondition
{
public:
    TwoFluidIsothermalSheathHagelaarBC(string name, string side) : BoundaryCondition(name, side){}
    ~TwoFluidIsothermalSheathHagelaarBC() {}
    virtual CellDataRef setBoundary();
    void getWallFlux(const CellDataRef innerU, CellDataRef ghostU);
};

#endif /* TwoFluidIsothermalSheathHagelaarBC_hpp */
