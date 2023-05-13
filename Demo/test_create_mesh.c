#include "hpc.h"


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    index brief = 0;

    mesh *refined_mesh;

    if (rank == 0) {
        index m = 2;
        index n = 3;

        mesh *newMesh = create_rect_mesh(m, n);
        //mesh_print(newMesh, brief);

        //printf("\nRefined mesh\n");

        refined_mesh = mesh_refine(newMesh);
        mesh_print(refined_mesh, brief);

        printf("\n");
    }

    mesh *local_mesh = scatter_mesh(refined_mesh, MPI_COMM_WORLD, rank, nof_processes);
    print_mpi(local_mesh, MPI_COMM_WORLD, rank, nof_processes, brief);

    MPI_Finalize();
}
