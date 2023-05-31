#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Skeleton {

	void Skeleton::Scatter(int rank, Mesh::LocalMesh &local_mesh) {

		// Broadcast Skeleton
		MPI_Bcast(&com_borders(0), (int) n_borders * 6,
                  MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&com_border_nodes.nodes(0), (int) com_border_nodes.nodes.count(),
                  MPI_LONG, 0, MPI_COMM_WORLD);
		
		CreateLocal(rank, local_mesh);
	}
}
#endif // _MPI
