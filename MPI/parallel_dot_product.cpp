#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Util {
    void parallel_dot_product(Vector<double> &v_acc, Vector<double> &v_dist, double &global_result) {
	assert(v_acc.count() == v_dist.count());
	// Compute local part of the dot product
	double local_result = 0;
	for (int i = 0; i < v_acc.count(); ++i) {
	    local_result += v_acc(i) * v_dist(i);
	}

	MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}
#endif // _MPI
