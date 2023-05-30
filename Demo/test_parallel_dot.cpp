#include <cstdio>
#include <iostream>
#include <vector>
#include <mpi.h>

#include "hpc.hpp"


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    std::vector<double> v_acc{1, 2, 3};
    std::vector<double> v_dist{(double) rank, (double) rank, (double) rank};
    
    double dot_result = 0;

    Util::parallel_dot_product(v_acc, v_dist, dot_result);

    if (rank == 0) {
	std::cout << "dot product result: " << dot_result << "\n";
    }

    MPI_Finalize();
}

