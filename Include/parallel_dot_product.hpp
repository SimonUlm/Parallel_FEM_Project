#ifndef HPC2_PARALLEL_DOT_PRODUCT_HPP
#define HPC2_PARALLEL_DOT_PRODUCT_HPP

#ifdef _MPI

#include <mpi.h>
#include <vector>

namespace Util {
    void ParallelDot(Vector < double > &v_acc, Vector < double > &v_dist, double & global_result);
}

#endif //_MPI

#endif // HPC2_PARALLEL_DOT_PRODUCT_HPP 
