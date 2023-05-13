#include "hpc.h"
#include <mpi.h>

void print_mpi(mesh *mesh_data, MPI_Comm comm, int rank, int nof_processes, int brief) {
    int dummy = 0;
    MPI_Status status;

    // Wait for previous process
    if (rank != 0) {
        MPI_Recv(&dummy, 1, MPI_INT, rank - 1, 0, comm, &status);
    }

    // Print data
    mesh_print(mesh_data, brief);

    // Signal next process
    if (rank != nof_processes - 1) {
        MPI_Send(&dummy, 1, MPI_INT, rank + 1, 0, comm);
    }

}