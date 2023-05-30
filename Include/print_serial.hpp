#ifndef HPC2_PRINT_SERIAL_HPP
#define HPC2_PRINT_SERIAL_HPP

#ifdef _MPI
#include <functional>
#include <mpi.h>

namespace MPI {
    void PrintSerial(MPI_Comm comm, int rank, int nof_processes, const std::function<void()> &print);
}
#endif

#endif //HPC2_PRINT_SERIAL_HPP
