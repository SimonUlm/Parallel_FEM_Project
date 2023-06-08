#ifndef HPC2_UTIL_MATRIX_HPP
#define HPC2_UTIL_MATRIX_HPP

#include <algorithm>
#include <cstdio>

#ifndef NDEBUG

#include <cassert>

#endif

#include "hpc.hpp"


namespace Util {

    enum StorageOrder {
        ROWMAJOR, COLMAJOR
    };

    /*
     * SedMatrix
     *
     * Matrix class for spares matrices in compress format with extracted diagonal
     *
     * The SedMatrix class is a storage format for sparse matrices. The matrix
     * can be only compressed column. It stores the diagonal in extracted form and
     * does not store zero entries
     *
     */
    class SedMatrix {
    public:
        /*
         * Constructor for nxn sparse matrix with nzmax non-zero elements
         *
         * n_: Number of columns/rows
         * nzmax_: Maximum number of nonzero elements
         *
         */
        SedMatrix(long n, long nzmax, bool is_symmetry_format = false) :
                n_(n), nzmax_(nzmax + 1),
                ptr_ind_(new long[nzmax + 1]()),
                data_(new double[nzmax + 1]()),
                is_symmetry_format_(is_symmetry_format) { ptr_ind_[n] = nzmax; }

        // Destructor
        ~SedMatrix() {
            delete[] ptr_ind_;
            delete[] data_;
        }

        SedMatrix(SedMatrix &&other) noexcept:
                ptr_ind_(other.ptr_ind_), data_(other.data_), nzmax_(other.nzmax_), n_(other.n_) {
            other.ptr_ind_ = nullptr;
            other.data_ = nullptr;
            other.nzmax_ = 0;
            other.n_ = 0;
        }

        SedMatrix(const SedMatrix &) = delete;

        // Assignment operations
        SedMatrix &operator=(SedMatrix &&other) = delete;

        SedMatrix &operator=(const SedMatrix &) = delete;

        // Data access operator with absolute index_ in array (direct access)
        const double &operator()(long i) const {
#ifndef NDEBUG
            assert(i < nzmax_);
#endif
            return data_[i];
        }

        double &operator()(long i) {
#ifndef NDEBUG
            assert(i < nzmax_);
#endif
            return data_[i];
        }

        // Data access operator with matrix coordinates (traversing through column need)
        double operator()(long i, long j) {
            // Check diagonal
            if (i == j) return data_[i];

            // Column
            long col_start = ptr_ind_[j];
            long col_end = ptr_ind_[j + 1];

            if (col_start == col_end) return 0;

            // Traverse col
            for (long k = col_start; k < col_end; ++k) {
                if (ptr_ind_[k] == i) {
                    return data_[k];
                }
            }

            return 0;
        }

        // Getter methods
        const long nzmax() const { return nzmax_; }

        const long n() const { return n_; }

        bool is_symmetry_format() { return is_symmetry_format_; }

        // Get ptr_ind from absolute position in array
        long get_ptr(long k) { return ptr_ind_[k]; }

        // Set ptr_ind at absolute position in array
        void set_ptr(long k, long ptr) { ptr_ind_[k] = ptr; }

        // Adding given value to matrix entry
        void add_val(long i, long j, double val);

        // Set column to zero except diagonal entry
        void zero_col(long j);

        // Set row to zero except diagonal entry
        void zero_row(long i);

        // Symmetric sed general matrix vector product
        // y <- alpha * A * x + beta * y
        void SymSpmv(double alpha, BlasVector &x, double beta, BlasVector &y);

        // Get copy of diagonal
        Util::BlasVector Diag() const {
            Util::BlasVector diag(n_);
            for (long i = 0; i < n_; ++i)
                diag(i) = data_[i];
            return diag;
        }

        void Init() {
            // Init Diagonal with 00 11 22 33 ..
            double v = 0;
            for (long k = 0; k < n_; ++k) {
                (*this)(k) = v;
                v += 11;
            }

            v = 0;
            long index = n_ + 1;
            ptr_ind_[0] = index;
            for (long k = 0; k < n_; ++k) {
                for (long j = 0; j < n_; ++j) {
                    if (j != k) {
                        data_[index] = v;
                        ptr_ind_[index] = j;

                        index += 1;
                    }
                    v += 10;
                }
                ptr_ind_[k + 1] = index;
                v = (double) k + 1;
            }

        }

        void InitDenseSpd() {
            long row_index = n_ + 1;

            // Set pointer and row indices for first column
            ptr_ind_[0] = row_index;
            for (long i = 1; i < n_ - 1; ++i)
                ptr_ind_[row_index++] = i;

            // Set all pointers and row indices
            for (long j = 1; j < n_ - 1; ++j) {
                ptr_ind_[j] = row_index;
                for (long i = 0; i < n_; ++i) {
                    if (i == j)
                        continue;
                    ptr_ind_[row_index++] = i;
                }
            }

            // Set pointer and row indices for last column
            ptr_ind_[n_ - 1] = row_index;
            for (long i = 1; i < n_ - 1; ++i)
                ptr_ind_[row_index++] = i;
            ptr_ind_[n_] = row_index;

            // Set values
            for (long i = 0; i < n_; ++i) {
                for (long j = 0; j < n_; ++j) {
                    auto value = (double) (n_ - std::abs(i - j) - 1);
                    if (value == 0)
                        continue;
                    add_val(i, j, value);
                }
            }
        }

        // Print data_ of matrix
        void Print();

    private:
        long nzmax_;                        // Maximum number of entries
        long n_;                            // Number of rows/columns
        long *ptr_ind_;                     // Col pointers and row indices
        double *data_;                      // Numerical values
        bool is_symmetry_format_ = false;   // If true, only half of the symmetric matrix is actually stored
    };

    /*
     * GeMatrix
     *
     * Matrix class for general dense matrices
     *
     * The GeMatrix class is a general storage format for dense matrices which stores
     * all values including zero values. The storage format can be row or column major
     */
    class GeMatrix {
    public:
        GeMatrix(GeMatrix &&) = delete;

        GeMatrix(const GeMatrix &) = delete;

        GeMatrix &operator=(GeMatrix &&other) = default;

        GeMatrix &operator=(const GeMatrix &) = delete;

        const double &operator()(long i, long j) const {
#ifndef NDEBUG
            assert(i < m_ && j < n_);
#endif
            return data_[i * inc_row_ + j * inc_col_];
        }

        double &operator()(long i, long j) {
#ifndef NDEBUG
            assert(i < m_ && j < n_);
#endif
            return data_[i * inc_row_ + j * inc_col_];
        }

        /*
         *   Constructor for GeMatrix
         *
         *   m_: Number of rows
         *   n_: Number of columns
         *   order: Storage order ROWMAJOR | COLMAJOR
         */
        GeMatrix(long m, long n, StorageOrder order) :
                m_(m), n_(n),
                inc_row_(order == StorageOrder::COLMAJOR ? 1 : n),
                inc_col_(order == StorageOrder::ROWMAJOR ? 1 : m),
                data_(new double[m * n]) {}

        /*
         *   Constructor for GeMatrix from SedMatrix (always as colmajor)
         *
         *   sed: Reference to SedMatrix
         *	 sym: Is given SedMatrix symmetric
         */
        explicit GeMatrix(SedMatrix &sed) :
                m_(sed.n()), n_(sed.n()),
                inc_row_(1), inc_col_(n_),
                data_(new double[sed.n() * sed.n()]) {
            FromSed(sed);
        }

        // Destructor
        ~GeMatrix() {
            delete[] data_;
        }

        // Print data_ of matrix
        void Print();

    private:
        long m_;       // Number of Rows
        long n_;       // Number of Columns
        long inc_row_;  // Row Stride
        long inc_col_;  // Column Stride
        double *data_; // Numerical values

        // Initialize Matrix from Sed Matrix
        void FromSed(SedMatrix &sed);
    };
}

#endif //HPC2_UTIL_MATRIX_HPP