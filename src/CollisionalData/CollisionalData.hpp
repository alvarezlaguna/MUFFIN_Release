//
//  CollisionalData.hpp
//  Muffin
// Just a singleton to access the mixture data
//
//  Created by Alejandro Alvarez Laguna on 17/02/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef CollisionalData_hpp
#define CollisionalData_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

class CollisionalData
{
public:
    static CollisionalData& getInstance()
    {static CollisionalData cd; return cd;}
    
    template <class T> void createData(const string& name)
    {data[name] = (void*)new T;}
    
    template <class T> void deleteData(const string& name)
    {delete (vector<T>*)data.find(name)->second;}
    
    template <class T> T& getData(const string& name)
    {return *(T*)data.find(name)->second;}
    
private:
    CollisionalData() {}
    CollisionalData(CollisionalData&) {}
    
    map<string, void*> data;
};

#endif /* CollisionalData_hpp */
