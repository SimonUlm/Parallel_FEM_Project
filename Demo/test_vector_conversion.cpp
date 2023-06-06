#include <iostream>
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
    int refine_factor = 1;

    GlobalMesh global_mesh;
    LocalMesh local_mesh;

    if (rank == 0)
    {
        global_mesh = GlobalMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine(refine_factor);
    }
    Skeleton::Skeleton skeleton(m, n, refine_factor, MPI_COMM_WORLD, rank);
    if (rank == 0)
        skeleton.Create(global_mesh);

    global_mesh.Scatter(local_mesh, skeleton);

    Util::Vector<double> accum_to_distr(local_mesh.get_n_nodes());
    accum_to_distr.Init(10);
    skeleton.AccumulatedToDistributed(accum_to_distr);

    Util::Vector<double> distr_to_accum(local_mesh.get_n_nodes());
    distr_to_accum.Init();
    skeleton.DistributedToAccumulated(distr_to_accum);

    MPI::PrintSerial(MPI_COMM_WORLD, rank, nof_processes, [&]() {
        std::cout << "Rank = " << rank << std::endl;
        for (auto &x : accum_to_distr)
            std::cout << x << " ";
        std::cout << std::endl;
        for (auto &x : distr_to_accum)
            std::cout << x << " ";
        std::cout << std::endl;
    });

    // Gather accum_to_distr onto root process
    Util::Vector<double> global_vector;
    Util::Vector<double> test_vector(local_mesh.get_n_nodes());
    test_vector.Init(50);
    skeleton.GatherAccumulatedVector(test_vector, global_vector);
    if (rank == 0)
        global_vector.Print();

    MPI_Finalize();
}
