#include <cstdio>
#include <stdexcept>

#include "hpc.hpp"

namespace Util {
	
	/*
	 * Adding given value to matrix entry
     *
     * i: Row index
     * j: Column index
     * val: value to be added
     *
	 */
	void SedMatrix::add_val(long i, long j, double val) {
		// Check if on diagonal
    	if (i == j) {
    		data[i] += val;
    		return;
    	 }
    	
    	// Column
    	long col_start = ptr_ind[j];
    	long col_end = ptr_ind[j+1];
		
		if (col_start == col_end) {
			throw std::invalid_argument(
				"Write position must not be zero entry! No memory was reserved for the given matrix entry"
			);
		}
		
		// Traverse column to find row
		for (long k = col_start; k < col_end; ++k) {
			if (ptr_ind[k] == i) {
				// row found
				data[k] += val;
				return;
			}
		}
		
		throw std::invalid_argument(
			"Write position must not be zero entry! No memory was reserved for the given matrix entry"
		);
	}
}