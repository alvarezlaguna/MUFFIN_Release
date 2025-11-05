//
//  TwoFluidIsothermalSheathAPHagelaarBC.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 29/01/19.
//  Copyright © 2019 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermalSheathAPHagelaarBC_hpp
#define TwoFluidIsothermalSheathAPHagelaarBC_hpp

#include <stdio.h>
#include "BoundaryCondition.hpp"

class TwoFluidIsothermalSheathAPHagelaarBC : public BoundaryCondition
{
public:
    TwoFluidIsothermalSheathAPHagelaarBC(string name, string side) : BoundaryCondition(name, side){}
    ~TwoFluidIsothermalSheathAPHagelaarBC() {}
    virtual CellDataRef setBoundary();
    void getWallFlux(const CellDataRef innerU, CellDataRef ghostU);
};

#endif /* TwoFluidIsothermalSheathAPHagelaarBC_hpp */
