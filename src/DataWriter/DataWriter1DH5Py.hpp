// DataWriter1DH5Py.hpp
// Writes 1D data by calling Python's h5py (via pybind11)

#ifndef DataWriter1DH5Py_hpp
#define DataWriter1DH5Py_hpp

#include "DataWriter1D.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class DataWriter1DH5Py : public DataWriter1D {
public:
    DataWriter1DH5Py(const std::string &name) : DataWriter1D(name) {}
    ~DataWriter1DH5Py() {}
    virtual void writeData(const int iter) override;
};

#endif // DataWriter1DH5Py_hpp
