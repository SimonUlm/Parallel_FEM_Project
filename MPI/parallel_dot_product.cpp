#include <cstdio>
#include <memory>
#include <vector>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Util {
    void parallel_dot_product(std::vector<long> &v_acc, std::vector<long> &v_dist, long &global_result) {
	// Compute local part of the dot product
	long local_result = 0;
	for (int i = 0; i < v_acc.size(); ++i) {
	    local_result += v_acc[i] * v_dist[i];
	}

	MPI_Allreduce(&local_result, global_result, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
}
#endif // _MPI
