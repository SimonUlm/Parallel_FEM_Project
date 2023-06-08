#include "hpc.hpp"

namespace Util {

    /*
     * Initialize Matrix from Sed Matrix
     *
     * sed: Reference to SedMatrix
     *
     */
    void GeMatrix::FromSed(SedMatrix &sed) {
        // Diagonal
        for (long i = 0; i < n_; ++i) {
            (*this)(i, i) = sed(i);
        }

        // Super-/Sub-Diagonal
        if (sed.is_symmetry_format()) {
            for (long j = 0; j < n_; ++j) {
                for (long p = sed.get_ptr(j); p < sed.get_ptr(j + 1); ++p) {
                    (*this)(sed.get_ptr(p), j) = sed(p);
                    (*this)(j, sed.get_ptr(p)) = sed(p);
                }
            }
        } else {
            for (long j = 0; j < n_; ++j) {
                for (long p = sed.get_ptr(j); p < sed.get_ptr(j + 1); ++p) {
                    (*this)(sed.get_ptr(p), j) = sed(p);
                }
            }
        }
    }
}
