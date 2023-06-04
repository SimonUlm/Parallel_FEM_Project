#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Skeleton {

	void Skeleton::Scatter(int rank, Mesh::LocalMesh &local_mesh) {

		// Broadcast Skeleton
		MPI_Bcast(com_borders.data(), (int) n_borders * 6,
                  MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(com_border_nodes.get_nodes().data(), (int) com_border_nodes.get_nodes().count(),
                  MPI_LONG, 0, MPI_COMM_WORLD);
		
		CreateLocal(rank, local_mesh);
	}
}
#endif // _MPI
