//
//  MeshData.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef MeshData_hpp
#define MeshData_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


using namespace std;
namespace py = pybind11;

class MeshData
{
public:
    static MeshData& getInstance()
    {static MeshData md; return md;}
    
    template <class T> void createData(const string& name, const int& size)
    {data[name] = (void*)new vector<T>(size);}

    template <class T> void create2DData(const string& name, const int& rows, const int& cols)
    {
        auto* arr = new py::array_t<T>({rows, cols}, {cols*sizeof(T), sizeof(T)});
        std::fill_n(arr->mutable_data(), rows * cols, T(0));
        data[name] = (void*)arr;
    }
    
    template <class T> void deleteData(const string& name)
    {delete (vector<T>*)data.find(name)->second;}
    
    template <class T> vector<T>& getData(const string& name)
    {return *(vector<T>*)data.find(name)->second;}
       
    template <class T> py::array_t<T>& get2DData(const string& name)
    {return *(py::array_t<T>*)data.find(name)->second;}

    template <class T> void delete2DData(const string& name)
    {delete (py::array_t<T>*)data.find(name)->second;}

    // Safe check whether a data entry exists
    bool hasData(const string &name) const
    {
        return data.find(name) != data.end();
    }


    
private:
    MeshData() {}
    MeshData(MeshData&) {}
    
    map<string, void*> data;
};

#endif /* MeshData_hpp */
