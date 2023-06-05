#ifndef HPC2_UTIL_MATRIX_HPP
#define HPC2_UTIL_MATRIX_HPP

#include <algorithm>
#include <cassert>
#include <cstdio>

#include "hpc.hpp"


namespace Util {	
	enum StorageOrder {ROWMAJOR, COLMAJOR};

	/* SedMatrix */
    
	class SedMatrix {
	/* 
     * Matrix class for spares matrixes in compress format with extraxted diagonal
     *
     * The SedMatrix class is a storage format for sparses matrixes which. The matrix
     * can be only compressed column. It stores the diagonal in extracted form and
     * does not store zero entries 
     * 
     */	
   	private:
   		long nzmax;			                // Maximum number of entries
   		long n;				                // Number of rows/columns
   		long *ptr_ind;		                // Col pointers and row indices
   		double *data;		                // Numerical values
        bool is_symmetry_format_ = false;   // If true, only half of the symmetric matrix is actually stored
   	public:
    	/*
    	 * Constructor for nxn sparse matrix with nzmax non-zero elements
         *
         * n: Number of columns/rows
         * nzmax: Maximum number of nonzero elements
         *
    	 */
    	SedMatrix(long n, long nzmax, bool is_symmetry_format = false) :
    			n(n), nzmax(nzmax+1), 
    			ptr_ind(new long[nzmax+1]()),
    			data(new double[nzmax+1]()),
                is_symmetry_format_(is_symmetry_format) {ptr_ind[n] = nzmax;}
    			
    	// Destructor
        ~SedMatrix() {
        	delete[] ptr_ind;
            delete[] data;
        }

        SedMatrix(SedMatrix &&other) noexcept:
                ptr_ind(other.ptr_ind), data(other.data), nzmax(other.nzmax), n(other.n) {
            other.ptr_ind = nullptr;
            other.data = nullptr;
            other.nzmax = 0;
            other.n = 0;
        }
        SedMatrix(const SedMatrix &) = delete;

        // Assignment operations
        SedMatrix &operator=(SedMatrix &&other) = delete;
        SedMatrix &operator=(const SedMatrix &) = delete;

        // Data access operator with absolute index in array (direct access)
        double & operator()(long i) const {
            assert(i < nzmax);
            return data[i];
        }
        double & operator()(long i) {
            assert(i < nzmax);
            return data[i];
        }

        // Data access operator with matrix coordinates (traversing through column need)
        double operator()(long i, long j) {
            // Check diagonal
            if (i == j) return data[i];

            // Column
            long col_start = ptr_ind[j];
            long col_end = ptr_ind[j+1];

            if (col_start == col_end) return 0;

            // Traverse col
            for (long k = col_start; k < col_end; ++k) {
                if (ptr_ind[k] == i) {
                    return data[k];
                }
            }

            return 0;
        }

    	// Getter methods		
    	long get_n(){return n;}
    	long get_nzmax(){return nzmax;}
        bool is_symmetry_format(){return is_symmetry_format_;}
    	
    	// Get ptr_ind from absolute position in array
    	long get_ptr(long k){return ptr_ind[k];}
    	
    	// Set ptr_ind at absolute postion in array
    	void set_ptr(long k, long ptr){ptr_ind[k] = ptr;}
    	
    	// Adding given value to matrix entry
    	void add_val(long i, long j, double val);
    	
    	// set column to zero except diagonal entry
    	void zero_col(long j);
    	
    	// set row to zero except diagonal entry
    	void zero_rows(long i, long j);
    	
    	// symmetric sed general matrix vector product
		// y <- alpha * A * x + beta * y 
        void SymSpmv(double alpha, BlasVector &x, double beta, BlasVector &y);

        // Get copy of diagonal
        Util::BlasVector Diag() const {
            Util::BlasVector diag(n);
            for (long i = 0; i < n; ++i)
                diag(i) = data[i];
            return diag;
        }
        
		void Init() {
    		// Init Diagonal with 00 11 22 33 ..
    		double v = 0;
    		for (long k = 0; k < n; ++k) {
    			(*this)(k) = v;
    			v += 11;	
    		}
    		
    		v = 0;
    		long index = n + 1;
    		ptr_ind[0] = index;
    		for (long k = 0; k < n; ++k) {
    			for (long j = 0; j < n; ++j) {
    				if (j != k) {
    					data[index] = v;
    					ptr_ind[index] = j;
    					
    					index += 1;
    				}
    				v += 10;
    			}
    			ptr_ind[k+1] = index;
    			v = k + 1;
    		}
    				
    	}

        void InitDenseSpd() {
            long row_index = n + 1;

            // Set pointer and row indices for first column
            ptr_ind[0] = row_index;
            for (long i = 1; i < n-1; ++i)
                ptr_ind[row_index++] = i;

            // Set all pointers and row indices
            for (long j = 1; j < n-1; ++j) {
                ptr_ind[j] = row_index;
                for (long i = 0; i < n; ++i) {
                    if (i == j)
                        continue;
                    ptr_ind[row_index++] = i;
                }
            }

            // Set pointer and row indices for last column
            ptr_ind[n-1] = row_index;
            for (long i = 1; i < n-1; ++i)
                ptr_ind[row_index++] = i;
            ptr_ind[n] = row_index;

            // Set values
            for (long i = 0; i < n; ++i) {
                for (long j = 0; j < n; ++j) {
                    auto value = (double) (n - std::abs(i - j) - 1);
                    if (value == 0)
                        continue;
                    add_val(i, j, value);
                }
            }
        }

        // Print data of matrix
    	void Print();
   	};
	
	/* GeMatrix */
    
    class GeMatrix {
    /*
     * Matrix class for general dense matrixes
     *
     * The GeMatrix class is a general storage format for dense matrixes which stores
     * all values including zero values. The storage format can be row or column major
     */	
    private:
    	long m;       // Number of Rows
    	long n;       // Number of Columns
    	long incRow;  // Row Stride
    	long incCol;  // Column Stride
        double *data; // Numerical values
        
        // Initialize Matrix from Sed Matrix
        void FromSed(SedMatrix &sed);

    public:
    	GeMatrix(GeMatrix &&) = delete;
        GeMatrix(const GeMatrix &) = delete;
        
        GeMatrix & operator=(GeMatrix &&other) = default;
        GeMatrix & operator=(const GeMatrix &) = delete;
    			
    	double & operator()(long i, long j) const {
    		assert(i < m && j < n);
            return data[i*incRow + j*incCol];
        }
        double & operator()(long i, long j) {
        	assert(i < m && j < n);
            return data[i*incRow + j*incCol];
        }
    
    	/*
    	 *   Constructor for GeMatrix
         *
         *   m: Number of rows
         *   n: Number of columns
         *   order: Storage order ROWMAJOR | COLMAJOR
    	 */
    	GeMatrix(long m, long n, StorageOrder order) :
    			m(m), n(n), 
    			incRow(order == StorageOrder::COLMAJOR? 1: n),
    			incCol(order == StorageOrder::ROWMAJOR? 1: m),
    			data(new double[m*n]) {}
    	
    	/*
    	 *   Constructor for GeMatrix from SedMatrix (always as colmajor)
         *
         *   sed: Reference to SedMatrix
         *	 sym: Is given SedMatrix symmetric
    	 */		
    	GeMatrix(SedMatrix &sed) :
    			m(sed.get_n()), n(sed.get_n()),
    			incRow(1), incCol(n),
    			data(new double[sed.get_n() * sed.get_n()]) {
            FromSed(sed);
		}
		
		// Destructor
        ~GeMatrix() {
            delete[] data;
        }
    	
		// Print data of matrix
        void Print();
   	};
}

#endif //HPC2_UTIL_MATRIX_HPP
