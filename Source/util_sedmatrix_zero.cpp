#include "hpc.hpp"

namespace Util {

    /*
	 * Set column j to zero except the diagonal entry
	 *
	 * j: column index
	 *
	 */
    void SedMatrix::zero_col(long j) {
        long col_start = ptr_ind_[j];
        long col_end = ptr_ind_[j + 1];

        // Traverse column to find row
        for (long k = col_start; k < col_end; ++k) {
            // row found
            data_[k] = 0;
        }
    }

    /*
     * Set the given row to zero except the diagonal entry
     *
     * i: row index
     *
     */
    void SedMatrix::zero_row(long i) {
        long col_start;
        long col_end;

        for (long k = 0; k < n_; ++k) {
            col_start = ptr_ind_[k];
            col_end = ptr_ind_[k + 1];

            // Traverse column to find row
            for (long l = col_start; l < col_end; ++l) {
                if (ptr_ind_[l] == i) {
                    // row found
                    data_[l] = 0;
                    break;
                }
            }
        }
    }
}
	
	
