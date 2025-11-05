//
//  LimiterFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LimiterFactory_hpp
#define LimiterFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>
#include "Limiter.hpp"

using namespace std;

// Factory for creating instances of Limiter (self-registration enabled)
class LimiterFactory
{
private:
    LimiterFactory();
    LimiterFactory(const LimiterFactory &) { }
    LimiterFactory &operator=(const LimiterFactory &) { return *this; }

    typedef std::function<std::unique_ptr<Limiter>(const std::string&)> CreateLimiterFn;
    std::map<std::string, CreateLimiterFn> m_registry;
    std::mutex m_registryMutex;
public:
    static LimiterFactory *Get()
    {
        static LimiterFactory instance;
        return &instance;
    }

    // Register a limiter creator. Returns true on success, false if name is already registered.
    static bool Register(const std::string &name, CreateLimiterFn pfnCreate);

    // Create limiter by name; returns a default limiter if name not found.
    static std::unique_ptr<Limiter> CreateLimiter(const std::string &name);
};

#endif /* LimiterFactory_hpp */
