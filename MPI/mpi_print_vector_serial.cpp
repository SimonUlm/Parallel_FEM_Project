#include <iostream>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Mesh;

void MpiPrintVectorSerial(std::vector<long> &vector_1, std::vector<long> &vector_2, MPI_Comm comm, int rank, int nof_processes) {
    int dummy = 0;
    MPI_Status status;

    // Wait for previous process
    if (rank != 0) {
        MPI_Recv(&dummy, 1, MPI_INT, rank - 1, 0, comm, &status);
    }

    // Print data
    std::cout << "Rank = " << rank << std::endl;
    for (long x : vector_1)
        std::cout << x << " ";
    std::cout << std::endl;
    for (long x : vector_2)
        std::cout << x << " ";
    std::cout << std::endl;

    // Signal next process
    if (rank != nof_processes - 1) {
        MPI_Send(&dummy, 1, MPI_INT, rank + 1, 0, comm);
    }
}
#endif // _MPI
