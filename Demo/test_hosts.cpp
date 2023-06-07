#include "hpc.hpp"
#include <unistd.h>
#include <mpi.h>
#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    printf("I am rank %d\n", rank);
    MPI_Finalize();
}
