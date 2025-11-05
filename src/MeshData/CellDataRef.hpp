//
//  CellDataRef.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//


#ifndef CellDataRef_hpp
#define CellDataRef_hpp

#include <stdio.h>
#include "../Parameters.hpp"
#include <pybind11/numpy.h>
namespace py = pybind11;

class CellDataRef
{
public:
    CellDataRef() : m_array(nullptr), m_iCell(0) {}
    CellDataRef(py::array_t<double>* array, int iCell, int offset = 0) : m_array(array), m_iCell(iCell), m_offset(offset) {};    
    CellDataRef withOffset(int offset) const {return CellDataRef(m_array, m_iCell,offset);}
    py::array_t<double> getView() const {
        py::buffer_info info = m_array->request();
        size_t itemSize = info.itemsize;
        size_t numRows = info.shape[0];
        size_t numCols = info.shape[1];
        size_t stride = info.strides[1];

        double* dataPtr = static_cast<double*>(info.ptr);
        double* columnPtr = dataPtr + m_iCell * stride / itemSize;

        return py::array_t<double>(
            {numRows},
            {numCols * itemSize},
            columnPtr
        );
    }
    double operator [](int i) const {return m_array->unchecked<2>()(i+m_offset, m_iCell);}
    double & operator [](int i) {return m_array->mutable_unchecked<2>()(i+m_offset, m_iCell);}
//     operator std::vector<double>() const {
//         std::vector<double> v(Parameters::NBEQS);
//         for (int iEq = 0; iEq < Parameters::NBEQS; iEq++)
//         {
//             v[iEq] = m_array->at(iEq, m_iCell);
//         }
//         return v;
// }
    
protected:
    py::array_t<double>* m_array;
    int m_iCell;
    int m_offset;
};

#endif /* CellDataRef_hpp */
