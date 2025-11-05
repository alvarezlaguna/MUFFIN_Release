//
//  BoundaryCondition.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 23/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef BoundaryCondition_hpp
#define BoundaryCondition_hpp

#include <stdio.h>
#include "../Component.hpp"
#include "../Parameters.hpp"
#include "../MeshData/CellDataRef.hpp"
#include "../MeshData/MeshData.hpp"
using namespace Parameters;

class BoundaryCondition : public Component
{
public:
    BoundaryCondition(string name, string side) : Component(name), m_side(side) {
        pybind11::array_t<double>& cells = MeshData::getInstance().get2DData<double>("Cells_py_CC");
        if (m_side == "Left")
        {
            m_value = CellDataRef(&cells, 0);
        }
        else if (m_side == "Right")
        {
            m_value = CellDataRef(&cells, NBCELLS+1);
        }
    }
    virtual ~BoundaryCondition(){}
    virtual CellDataRef setBoundary() = 0;
protected:
    CellDataRef m_value;
    string m_side;
};

#endif /* BoundaryCondition_hpp */
