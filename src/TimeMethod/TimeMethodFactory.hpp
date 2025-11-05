//
//  TimeMethodFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TimeMethodFactory_hpp
#define TimeMethodFactory_hpp

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>
#include "TimeMethod.hpp"

using namespace std;

// Factory for creating instances of TimeMethod (self-registration enabled)
class TimeMethodFactory
{
private:
    TimeMethodFactory();
    TimeMethodFactory(const TimeMethodFactory &) { }
    TimeMethodFactory &operator=(const TimeMethodFactory &) { return *this; }

    typedef std::function<std::unique_ptr<TimeMethod>(const std::string&)> CreateTimeMethodFn;
    std::map<std::string, CreateTimeMethodFn> m_registry;
    std::mutex m_registryMutex;
public:
    static TimeMethodFactory *Get()
    {
        static TimeMethodFactory instance;
        return &instance;
    }

    // Register a time-method creator. Returns true on success, false if name already registered.
    static bool Register(const std::string &name, CreateTimeMethodFn pfnCreate);

    // Create a time-method by name; returns a default if name not found.
    static std::unique_ptr<TimeMethod> CreateTimeMethod(const std::string &name);
};

#endif /* TimeMethodFactory_hpp */
