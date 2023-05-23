#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Mesh;

void GlobalMesh::Scatter(LocalMesh &local_mesh, MPI_Comm comm, int rank, int nof_processes) {

    // Allocate global mesh on local level
    long mesh_data[5];
    if (rank == 0) {
        mesh_data[0] = m;
        mesh_data[1] = n;
        mesh_data[2] = nodes.count;
        mesh_data[3] = elements.count;
        mesh_data[4] = boundary.count;
    }
    MPI_Bcast(mesh_data, 5, MPI_LONG, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        m = mesh_data[0];
        n = mesh_data[1];
        nodes = List<Node>(mesh_data[2]);
        elements = List<Element>(mesh_data[3]);
        boundary = List<BoundaryEdge>(mesh_data[4]);
    }

    // Send global mesh to all nodes
    MPI_Bcast(&nodes(0).x, (int) nodes.count * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&elements(0).n1, (int) elements.count * 7, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boundary(0).n1, (int) boundary.count * 4, MPI_LONG, 0, MPI_COMM_WORLD);
}
#endif // _MPI
