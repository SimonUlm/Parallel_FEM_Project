#include <cstdio>
#include <mpi.h>

#include "mesh.hpp"

using namespace Mesh;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    int m = 3;
    int n = 2;

    RectangularMesh global_mesh;

    int nof_local_elem;
    if (rank == 0)
    {
        global_mesh = RectangularMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine();
        nof_local_elem = (int) global_mesh.elements.count / nof_processes;
    }
    MPI_Bcast(&nof_local_elem, 1, MPI_INT, 0, MPI_COMM_WORLD);

    RectangularMesh local_mesh(1, 1, 0, nof_local_elem, 0);
    global_mesh.Scatter(local_mesh, MPI_COMM_WORLD, rank, nof_local_elem);

    MpiPrintSerial(local_mesh, MPI_COMM_WORLD, rank, nof_processes);

    MPI_Finalize();
}
