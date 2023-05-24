#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Mesh;

void MpiPrintSerial(LocalMesh &mesh, Skeleton::Skeleton &skeleton, MPI_Comm comm, int rank, int nof_processes) {
    int dummy = 0;
    MPI_Status status;

    // Wait for previous process
    if (rank != 0) {
        MPI_Recv(&dummy, 1, MPI_INT, rank - 1, 0, comm, &status);
    }

    // Print data
    mesh.Print();
    skeleton.Print();

    // Signal next process
    if (rank != nof_processes - 1) {
        MPI_Send(&dummy, 1, MPI_INT, rank + 1, 0, comm);
    }
}
#endif // _MPI
