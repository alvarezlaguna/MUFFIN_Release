#ifndef SourceTermRegistrar_hpp
#define SourceTermRegistrar_hpp

#include <string>
#include <memory>
#include "SourceTermFactory.hpp"

// Helper to register a SourceTerm type T with the factory at static initialization time.
template<typename T>
struct SourceTermRegistrar {
    SourceTermRegistrar(const std::string &name) {
        SourceTermFactory::Register(name, [](const std::string &n)->std::unique_ptr<SourceTerm> {
            return std::unique_ptr<SourceTerm>(new T(n));
        });
    }
};

// Macro for convenient registration in a .cpp file
#define REGISTER_SOURCETERM(NAME, CLASS) \
    static SourceTermRegistrar<CLASS> registrar_##CLASS(NAME)

#endif /* SourceTermRegistrar_hpp */
