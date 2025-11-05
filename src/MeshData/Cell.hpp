//
//  Cell.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Cell_hpp
#define Cell_hpp

#include <stdio.h>

class Cell
{
public:
    Cell(const int ID){m_cellID = ID;}
    Cell(){}
    virtual ~Cell(){}
    
protected:
    int m_cellID;
};

#endif /* Cell_hpp */
