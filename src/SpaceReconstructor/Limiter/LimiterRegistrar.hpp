#ifndef LimiterRegistrar_hpp
#define LimiterRegistrar_hpp

#include <string>
#include <memory>
#include "LimiterFactory.hpp"

// Helper to register a limiter type T with the factory at static initialization time.
template<typename T>
struct LimiterRegistrar {
    LimiterRegistrar(const std::string &name) {
        LimiterFactory::Register(name, [](const std::string &n)->std::unique_ptr<Limiter> {
            return std::unique_ptr<Limiter>(new T(n));
        });
    }
};

// Macro for convenient registration in a .cpp file
#define REGISTER_LIMITER(NAME, CLASS) \
    static LimiterRegistrar<CLASS> registrar_##CLASS(NAME)

#endif /* LimiterRegistrar_hpp */
