

#ifndef InputFilesProvider_hpp
#define InputFilesProvider_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>


namespace py=pybind11;


class InputFilesProvider
{
public:
    static InputFilesProvider &getInstance(){
        static InputFilesProvider instance;
        return instance;
    };

    //InputFilesProvider(InputFilesProvider &other) = delete;
    void operator=(const InputFilesProvider &) = delete;

    void setData(std::string name,py::array_t<double> arr);


    py::array_t<double>& getData(std::string name);


protected:
    InputFilesProvider(/* args */);
    std::map<std::string, py::array_t<double>> data;


    

};

#endif
