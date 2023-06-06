#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Solver {
    double ParallelDot(Util::Vector<double> &v_acc, Util::Vector<double> &v_dist) {
	assert(v_acc.count() == v_dist.count());
	// Compute local part of the dot product
	double local_result = 0;
    double global_result = 0;
	for (int i = 0; i < v_acc.count(); ++i)
	    local_result += v_acc(i) * v_dist(i);

	MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_result;
    }
}
#endif // _MPI
