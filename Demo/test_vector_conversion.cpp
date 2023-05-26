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

    if (rank == 0)
    {
        global_mesh = GlobalMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine();
        global_mesh.Refine();
    }

    global_mesh.Scatter(local_mesh, MPI_COMM_WORLD, rank, nof_processes);

    std::vector<long> accum_to_distr(local_mesh.get_n_nodes(), 10);
    local_mesh.get_vector_converter().AccumulatedToDistributed(accum_to_distr);

    std::vector<long> distr_to_accum(local_mesh.get_n_nodes(), 10);
    local_mesh.get_vector_converter().DistributedToAccumulated(distr_to_accum, MPI_COMM_WORLD);

    MpiPrintVectorSerial(accum_to_distr, distr_to_accum, MPI_COMM_WORLD, rank, nof_processes);

    MPI_Finalize();
}
