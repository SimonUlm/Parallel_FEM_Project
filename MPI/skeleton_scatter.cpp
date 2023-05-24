#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Skeleton {

	void Skeleton::Scatter(int rank, Mesh::LocalMesh &local_mesh) {
		long n_nodes = comBorderNodes.get_n_nodes();
		long n_borders = get_n_borders();
		
		if (rank != 0)
			*this = Skeleton(n_borders, n_nodes, LOCAL);
		
		// Broadcast Skeleton
		MPI_Bcast(&comBorders(0), (int) n_borders * 6,
			  MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&comBorderNodes.nodes(0), (int) n_borders * n_nodes,
			  MPI_LONG, 0, MPI_COMM_WORLD);
		
		CreateLocal(rank, local_mesh);	
	}
}
#endif // _MPI
