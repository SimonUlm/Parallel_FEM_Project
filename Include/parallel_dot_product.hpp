#ifndef HPC2_PARALLEL_DOT_PRODUCT_HPP
#define HPC2_PARALLEL_DOT_PRODUCT_HPP

#include <mpi.h>
#include <vector>

namespace Util {
    void parallel_dot_product(Vector<double> &v_acc, Vector<double> &v_dist, double &global_result);
}

#endif // HPC2_PARALLEL_DOT_PRODUCT_HPP 
