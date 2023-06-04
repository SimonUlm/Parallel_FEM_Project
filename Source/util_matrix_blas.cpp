#include "hpc.hpp"

namespace Util {
    
    // ----- Level 2 -----
    
    // y <- alpha * A * x + beta * y
    void SedMatrix::SymSpmv(double alpha, BlasVector &x, double beta, BlasVector &y) {
        
		y.Scal(beta);
		
		if (alpha == 0)
            return;

        if (is_symmetry_format_) {
            // Traverse columns
            for (long j = 0; j < n; ++j) {
                // Diagonal
                y(j) += alpha * data[j] * x(j);

                // Sub/Super-Diagonal
                for (long p = ptr_ind[j]; p < ptr_ind[j + 1]; ++p) {
                    // lower part axpy-based
                    y(ptr_ind[p]) += alpha * data[p] * x(j);

                    // upper part ddot-based
                    y(j) += alpha * data[p] * x(ptr_ind[p]);
                }
            }
        } else {
            // Traverse columns
            for (long j = 0; j < n; ++j) {
                // Diagonal
                y(j) += alpha * data[j] * x(j);

                // Sub/Super-Diagonal
                for (long p = ptr_ind[j]; p < ptr_ind[j + 1]; ++p) {
                    // lower part axpy-based
                    y(ptr_ind[p]) += alpha * data[p] * x(j);
                }
            }
        }
	}		

}
