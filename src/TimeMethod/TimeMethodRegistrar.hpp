#ifndef TimeMethodRegistrar_hpp
#define TimeMethodRegistrar_hpp

#include <string>
#include <memory>
#include "TimeMethodFactory.hpp"

// Helper to register a TimeMethod type T with the factory at static initialization time.
template<typename T>
struct TimeMethodRegistrar {
    TimeMethodRegistrar(const std::string &name) {
        TimeMethodFactory::Register(name, [](const std::string &n)->std::unique_ptr<TimeMethod> {
            return std::unique_ptr<TimeMethod>(new T(n));
        });
    }
};

// Macro for convenient registration in a .cpp file
#define REGISTER_TIMEMETHOD(NAME, CLASS) \
    static TimeMethodRegistrar<CLASS> registrar_##CLASS(NAME)

#endif /* TimeMethodRegistrar_hpp */
