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
            for (long j = 0; j < n_; ++j) {
                // Diagonal
                y(j) += alpha * data_[j] * x(j);

                // Sub/Super-Diagonal
                for (long p = ptr_ind_[j]; p < ptr_ind_[j + 1]; ++p) {
                    // lower part axpy-based
                    y(ptr_ind_[p]) += alpha * data_[p] * x(j);

                    // upper part ddot-based
                    y(j) += alpha * data_[p] * x(ptr_ind_[p]);
                }
            }
        } else {
            // Traverse columns
            for (long j = 0; j < n_; ++j) {
                // Diagonal
                y(j) += alpha * data_[j] * x(j);

                // Sub/Super-Diagonal
                for (long p = ptr_ind_[j]; p < ptr_ind_[j + 1]; ++p) {
                    // lower part axpy-based
                    y(ptr_ind_[p]) += alpha * data_[p] * x(j);
                }
            }
        }
    }

}
