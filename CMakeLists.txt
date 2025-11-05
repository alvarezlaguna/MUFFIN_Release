cmake_minimum_required(VERSION 2.8.12)
project(muffin)

add_subdirectory(pybind11)

set(CPPFLAGS "-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_MODULE_PATH}/cmake")

find_package (Python3 COMPONENTS Interpreter Development REQUIRED)
include(FindPkgConfig)
pkg_check_modules(PETSC PETSc REQUIRED)

file(GLOB_RECURSE Muffin_SOURCES "src/*.cpp")
file(GLOB_RECURSE Muffin_HEADERS "src/*.h")

set (Muffin_INCLUDE_DIRS "")

foreach (_headerFile ${Muffin_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND Muffin_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES Muffin_INCLUDE_DIRS)

find_package(cxxopts REQUIRED)
message("Python3_SITELIB: ${Python3_SITELIB}")


pybind11_add_module(muffin ${Muffin_SOURCES})
target_link_libraries(muffin PRIVATE ${PETSC_LIBRARIES} ${MPI_LIBRARIES})
target_include_directories(muffin PRIVATE ${Muffin_INCLUDE_DIRS} )

# Install the library
install(TARGETS muffin
        LIBRARY DESTINATION ${Python3_SITELIB}
        ARCHIVE DESTINATION ${Python3_SITELIB}
        RUNTIME DESTINATION ${Python3_SITELIB})