//
//  Cell.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 14/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Cell1D_hpp
#define Cell1D_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "../Parameters.hpp"
#include "Cell.hpp"
#include "CellDataRef.hpp"

using namespace std;
using namespace Parameters;

class Cell1D : public Cell
{
public:
    Cell1D() : uCC(), uL(), uR() {}
    Cell1D(CellDataRef uCCRef, CellDataRef uLRef, CellDataRef uRRef) : uCC(uCCRef), uL(uLRef), uR(uRRef){};
    Cell1D(const int ID, CellDataRef uCCRef, CellDataRef uLRef, CellDataRef uRRef) : uCC(uCCRef), uL(uLRef), uR(uRRef) { m_cellID = ID; };

    void setCellID(const int ID) { m_cellID = ID; }

    void setCCValue(const double value) { uCC[0] = value; }
    void setCCValue(const vector<double> value)
    {
        for (int iEq = 0; iEq < Parameters::NBEQS; iEq++)
        {
            uCC[iEq] = value[iEq];
        }
    }

    void setLeftValue(const double value) { uL[0] = value; }
    void setLeftValue(const vector<double> value)
    {
        for (int iEq = 0; iEq < Parameters::NBEQS; iEq++)
        {
            uL[iEq] = value[iEq];
        }
    }

    void setRightValue(const double value) { uR[0] = value; }
    void setRightValue(const vector<double> value)
    {
        for (int iEq = 0; iEq < Parameters::NBEQS; iEq++)
        {
            uR[iEq] = value[iEq];
        }
    }

    CellDataRef u_CC() { return uCC; } // get function
    CellDataRef u_L() { return uL; }
    CellDataRef u_R() { return uR; }
    int ID() { return cellID; }

    ~Cell1D() {}

public:
    CellDataRef uCC; // Cell center value of the cell
    CellDataRef uR;  // Right value of the cell
    CellDataRef uL;  // Left value of the cell
    int cellID;      // ID of the cell
    double dx;
};

#endif /* Cell1D_hpp */
