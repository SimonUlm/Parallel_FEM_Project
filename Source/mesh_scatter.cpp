#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Mesh;

void  RectangularMesh::Scatter(RectangularMesh &local_mesh, MPI_Comm comm, int rank, int nof_local_elem) {

    if (rank == 0) {
        MPI_Scatter(&elements(0), nof_local_elem * 7, MPI_LONG,
                    &local_mesh.elements(0), nof_local_elem * 7, MPI_LONG, 0, comm);
    } else {
        MPI_Scatter(nullptr, 0, nullptr,
                    &local_mesh.elements(0), nof_local_elem * 7, MPI_LONG, 0, comm);
    }
}
#endif // _MPI
