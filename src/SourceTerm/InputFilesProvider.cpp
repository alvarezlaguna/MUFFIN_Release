






#include <stdio.h>
#include <map>
#include <string>
#include <vector>
#include "InputFilesProvider.hpp"
#include <pybind11/pybind11.h>

using namespace py;


InputFilesProvider::InputFilesProvider(/* args */)
{

}




void InputFilesProvider::setData(std::string name, py::array_t<double> arr){
    data.insert(std::pair<std::string, py::array_t<double>>(name,arr));

}


py::array_t<double>& InputFilesProvider::getData(std::string name){

    return (data.find(name)->second);
}

