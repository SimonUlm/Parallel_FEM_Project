#include <cstdio>
#include <memory>
#include <vector>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Util {
    void parallel_dot_product(std::vector<double> &v_acc, std::vector<double> &v_dist, double &global_result) {
	assert(v_acc.size() == v_dist.size());
	// Compute local part of the dot product
	double local_result = 0;
	for (int i = 0; i < v_acc.size(); ++i) {
	    local_result += v_acc[i] * v_dist[i];
	}

	MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}
#endif // _MPI
