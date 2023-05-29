#ifndef HPC2_PARALLEL_DOT_PRODUCT_HPP
#define HPC2_PARALLEL_DOT_PRODUCT_HPP

#include <mpi.h>
#include <vector>

void parallel_dot_product(std::vector<long>& v_acc, std::vector<long> &v_dist, long &global_result);

#endif // HPC2_PARALLEL_DOT_PRODUCT_HPP 
