#include <mpi.h>

#include "hpc.hpp"

namespace Skeleton {

    void Skeleton::Scatter(Mesh::LocalMesh &local_mesh, long m, long n, long refine_factor, MPI_Comm comm, int rank) {

        // If process is not root, skeleton is empty so far
        if (rank != 0)
            *this = Skeleton(m, n, refine_factor);

        // Broadcast Skeleton
        MPI_Bcast(com_borders_.data(), (int) n_borders_ * 6,
                  MPI_LONG, 0, comm_);
        MPI_Bcast(com_border_nodes_.nodes().data(), (int) com_border_nodes_.nodes().count(),
                  MPI_LONG, 0, comm_);

        CreateLocal(local_mesh, comm, rank);
    }
}

