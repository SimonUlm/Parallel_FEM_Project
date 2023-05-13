#include "hpc.h"
#include <mpi.h>

mesh * scatter_mesh(mesh *global_mesh, MPI_Comm comm, int rank, int nof_processes) {
    mesh *local_mesh;

    index a = 0;
    int b = a;

    // Broadcast basic information needed for local memory allocation
    int nof_local_elem;
    if (rank == 0)
        nof_local_elem = global_mesh->nelem / nof_processes;
    MPI_Bcast(&nof_local_elem, 1, MPI_INT, 0, comm);

    // Allocate mesh data structs
    local_mesh = mesh_alloc_with_edges(0, nof_local_elem, 0, 0, 0);

    if (rank == 0) {
        MPI_Scatter(global_mesh->elem, nof_local_elem * 7, MPI_INT,
                    local_mesh->elem, nof_local_elem * 7, MPI_INT, 0, comm);
    } else {
        MPI_Scatter(0, 0, 0,
                    local_mesh->elem, nof_local_elem * 7, MPI_INT, 0, comm);
    }

    return local_mesh;
}