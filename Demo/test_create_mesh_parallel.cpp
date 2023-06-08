#include <cstdio>
#include <mpi.h>

#include "hpc.hpp"

using namespace Mesh;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    int m = 3;
    int n = 2;

    GlobalMesh global_mesh;
    LocalMesh local_mesh;
    Skeleton::Skeleton skeleton;

    if (rank == 0)
    {
        global_mesh.Create(m, n);
        global_mesh.Refine();
        global_mesh.Refine();
        skeleton.Create(global_mesh, m, n);
    }

    global_mesh.Scatter(local_mesh, skeleton, MPI_COMM_WORLD, rank);

    MPI::PrintSerial(MPI_COMM_WORLD, rank, nof_processes, [&]() {
        local_mesh.Print();
        skeleton.Print();
    });

    MPI_Finalize();
}
