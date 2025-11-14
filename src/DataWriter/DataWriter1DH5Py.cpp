// DataWriter1DH5Py.cpp
// Implementation that calls Python (h5py + numpy) via pybind11 to write HDF5 files.

#include "DataWriter1DH5Py.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../Parameters.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <iostream>

using namespace Parameters;

void DataWriter1DH5Py::writeData(const int iter)
{
    // Gather data from MeshData
    std::vector<Cell1D>& cells  = MeshData::getInstance().getData<Cell1D>("Cells");
    std::vector<double>& x = MeshData::getInstance().getData<double>("x");
    double physTime = MeshData::getInstance().getData<double>("physTime")[0];

    m_iter = iter;

    char fname_c[512];
    snprintf(fname_c, sizeof(fname_c), "%s/file_iter_%06i_time_%.4e.h5", RESULTDIRECTORY.c_str(), m_iter, physTime);
    std::string fname(fname_c);

    // Convert data to numpy arrays
    py::gil_scoped_acquire acquire;
    try {
        py::module_ h5py = py::module_::import("h5py");
        py::module_ np = py::module_::import("numpy");

        // create numpy array for x (explicit buffer request)
        py::array_t<double> x_arr((size_t)x.size());
        {
            auto buf = x_arr.request();
            double *ptr = static_cast<double*>(buf.ptr);
            std::memcpy(ptr, x.data(), x.size()*sizeof(double));
        }

        // create 2D numpy array for u (NBCELLS x NBEQS)
        std::vector<double> buf(NBCELLS * NBEQS);
        for (unsigned int i = 0; i < NBCELLS; ++i){
            for (unsigned int j = 0; j < NBEQS; ++j){
                buf[i*NBEQS + j] = cells[i].u_CC()[j];
            }
        }
        py::array_t<double> u_arr({(size_t)NBCELLS, (size_t)NBEQS});
        {
            auto bufinfo = u_arr.request();
            double *ptr = static_cast<double*>(bufinfo.ptr);
            std::memcpy(ptr, buf.data(), buf.size()*sizeof(double));
        }

    // open file and write (use h5py assignment to avoid dimensionality issues)
    py::object f = h5py.attr("File")(fname, "w");
    f["x"] = x_arr;
    f["u"] = u_arr;

        // If Phi exists in MeshData, write it as well
        if (MeshData::getInstance().hasData("Phi")){
            try {
                std::vector<double>& phi = MeshData::getInstance().getData<double>("Phi");
                if (phi.size() == x.size()){
                    py::array_t<double> phi_arr((size_t)phi.size());
                    auto pbuf = phi_arr.request();
                    std::memcpy(pbuf.ptr, phi.data(), phi.size()*sizeof(double));
                    f["Phi"] = phi_arr;
                } else {
                    // size mismatch: still write if possible by using min size
                    size_t n = std::min(phi.size(), x.size());
                    py::array_t<double> phi_arr((size_t)n);
                    auto pbuf = phi_arr.request();
                    std::memcpy(pbuf.ptr, phi.data(), n*sizeof(double));
                    f["Phi"] = phi_arr;
                }
            }
            catch (const std::exception &e){
                std::cerr << "Error retrieving Phi from MeshData: " << e.what() << std::endl;
            }
        }

    // write time as dataset
    f["time"] = py::cast(physTime);
        f.attr("close")();
    }
    catch (py::error_already_set &e){
        std::cerr << "Python error when writing HDF5: " << e.what() << std::endl;
    }
}
