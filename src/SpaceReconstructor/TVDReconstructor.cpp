//
//  TVDReconstructor.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 13/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "TVDReconstructor.hpp"
#include "../Parameters.hpp"
#include <pybind11/pybind11.h>


namespace py = pybind11;
using namespace std;
using namespace Parameters;

void TVDReconstructor::setup() {
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // Initialize the pointer to the limiter with the name from the options
    m_limiter = LimiterFactory::CreateLimiter(LIMITERNAME);
    if(my_rank == MPI_WRITER)
        py::print("Reconstructor using limiter:\t ",m_limiter->getName(),"\n");

}

void TVDReconstructor::unsetup() {
    //m_limiter.release();
}
