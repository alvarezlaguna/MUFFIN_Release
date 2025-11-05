#ifndef ReconstructorRegistrar_hpp
#define ReconstructorRegistrar_hpp

#include <string>
#include <memory>
#include "ReconstructorFactory.hpp"

// Helper to register a reconstructor type T with the factory at static initialization time.
template<typename T>
struct ReconstructorRegistrar {
    ReconstructorRegistrar(const std::string &name) {
        ReconstructorFactory::Register(name, [](const std::string &n)->std::unique_ptr<SpaceReconstructor> {
            return std::unique_ptr<SpaceReconstructor>(new T(n));
        });
    }
};

// Macro for convenient registration in a .cpp file
#define REGISTER_RECONSTRUCTOR(NAME, CLASS) \
    static ReconstructorRegistrar<CLASS> registrar_##CLASS(NAME)

#endif /* ReconstructorRegistrar_hpp */
