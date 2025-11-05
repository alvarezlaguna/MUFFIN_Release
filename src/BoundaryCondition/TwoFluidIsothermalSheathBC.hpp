//
//  TwoFluidIsothermalSheathBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 03/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermalSheathBC_hpp
#define TwoFluidIsothermalSheathBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class TwoFluidIsothermalSheathBC : public BoundaryCondition
{
public:
    TwoFluidIsothermalSheathBC(string name, string side) : BoundaryCondition(name, side){}
    ~TwoFluidIsothermalSheathBC() {}
    virtual CellDataRef setBoundary();
    void getWallFlux(const CellDataRef innerU, CellDataRef ghostU);
};

#endif /* TwoFluidIsothermalSheathBC_hpp */
