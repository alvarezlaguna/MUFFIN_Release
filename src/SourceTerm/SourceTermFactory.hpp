//
//  SourceTermFactory.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef SourceTermFactory_hpp
#define SourceTermFactory_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include "SourceTerm.hpp"
#include <stdio.h>
#include <iostream>
#include <map>
#include <functional>
#include <memory>
#include <mutex>

using namespace std;

// Factory for creating instances of SourceTerm (self-registration enabled)
class SourceTermFactory
{
private:
    SourceTermFactory();
    SourceTermFactory(const SourceTermFactory &) { }
    SourceTermFactory &operator=(const SourceTermFactory &) { return *this; }

    typedef std::function<std::unique_ptr<SourceTerm>(const std::string&)> CreateSourceTermFn;
    std::map<std::string, CreateSourceTermFn> m_registry;
    std::mutex m_registryMutex;
public:
    static SourceTermFactory *Get()
    {
        static SourceTermFactory instance;
        return &instance;
    }

    // Register a source-term creator. Returns true on success, false if name already registered.
    static bool Register(const std::string &name, CreateSourceTermFn pfnCreate);

    // Create a source-term by name; returns a default if not found.
    static std::unique_ptr<SourceTerm> CreateSourceTerm(const std::string &name);
};

#endif /* SourceTermFactory_hpp */
