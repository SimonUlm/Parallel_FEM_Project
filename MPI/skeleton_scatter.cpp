#include <array>
#include <cstdio>
#include <memory>
#include <mpi.h>

#include "hpc.hpp"

namespace Skeleton {

	void Skeleton::Scatter(Mesh::LocalMesh &local_mesh) {

		// Broadcast Skeleton
		MPI_Bcast(com_borders.data(), (int) n_borders * 6,
                  MPI_LONG, 0, comm);
		MPI_Bcast(com_border_nodes.get_nodes().data(), (int) com_border_nodes.get_nodes().count(),
                  MPI_LONG, 0, comm);
		
		CreateLocal(local_mesh);
	}
} // Namespace Skeleton

