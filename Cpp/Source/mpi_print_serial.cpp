#include "mesh.hpp"
#include <mpi.h>

using namespace Mesh;

void mpi_print_serial(RectangularMesh &mesh, MPI_Comm comm, int rank, int nof_processes) {
    int dummy = 0;
    MPI_Status status;

    // Wait for previous process
    if (rank != 0) {
        MPI_Recv(&dummy, 1, MPI_INT, rank - 1, 0, comm, &status);
    }

    // Print data
    mesh.print();

    // Signal next process
    if (rank != nof_processes - 1) {
        MPI_Send(&dummy, 1, MPI_INT, rank + 1, 0, comm);
    }
}