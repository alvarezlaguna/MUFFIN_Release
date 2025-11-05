//
//  Options.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 01/03/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Options_hpp
#define Options_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

class Options
{
public:
    static Options& getInstance()
    {static Options cd; return cd;}
    
    template <class T> void createData(const string& name)
    {data[name] = (void*)new T;}
    
    template <class T> void deleteData(const string& name)
    {delete (vector<T>*)data.find(name)->second;}
    
    template <class T> T& getData(const string& name)
    {return *(T*)data.find(name)->second;}
    
private:
    Options() {}
    Options(Options&) {}
    
    map<string, void*> data;
};

#endif /* Options_hpp */
