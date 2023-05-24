#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Skeleton {

	void Skeleton::Scatter(int rank, Mesh::LocalMesh &local_mesh) {

		// Broadcast Skeleton
		MPI_Bcast(&comBorders(0), (int) n_borders * 6,
			  MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&comBorderNodes.nodes(0), (int) comBorderNodes.nodes.count,
			  MPI_LONG, 0, MPI_COMM_WORLD);
		
		CreateLocal(rank, local_mesh);
	}
}
#endif // _MPI
