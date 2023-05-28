#ifndef HPC2_UTIL_MATRIX_HPP
#define HPC2_UTIL_MATRIX_HPP

#include <cassert>

#include "hpc.hpp"


namespace Util {	
	enum StorageOrder {ROWMAJOR, COLMAJOR};
	
	/* SedMatrix */
    /*!  \class SedMatrix util_matrix.hpp "Include/util_matrix.hpp"
     *   \brief Matrix class for spares matrixes in compress format with extraxted diagonal
     *
     *	 The SedMatrix class is a storage format for sparses matrixes which. The matrix
     *	 can be either compressed row or column. It stores the diagonal in extracted form and
     *   does not store zero entries  
     */	
	class SedMatrix {
   	private:
   		long nzmax;			/*!< Maximum number of entries */
   		long n;				/*!< Number of rows/columns */
   		long *i;			/*!< Col pointers and row indices */
   		double *data;		/*!< Numerical values */	
   	public:
   		SedMatrix(SedMatrix &&) = delete;
        SedMatrix(const SedMatrix &) = delete;
        
        SedMatrix & operator=(SedMatrix &&other) = default;
        SedMatrix & operator=(const SedMatrix &) = delete;
        
        double & operator()(long i) const {
    		assert(i < nzmax);
            return data[i];
        }
        double & operator()(long i) {
        	assert(i < nzmax);
            return data[i];
        }
    	
    	SedMatrix(long n, long nzmax) :
    			n(n), nzmax(nzmax), 
    			i(new long[nzmax]),
    			data(new double[nzmax]) {i[nzmax-1] = nzmax-1;}
    			
    	long get_n(){return n;}
    	
    	long get_i(long k){return i[k];}
    	
    	void Print();
    	
    	// Debug purpose
    	void Init() {
    		// Init Diagonal with 00 11 22 33 ..
    		double v = 0;
    		for (long k = 0; k < n; ++k) {
    			(*this)(k) = v;
    			v += 11;	
    		}
    				
    	}
    	
    	// Missing ? sed *sed_nz_pattern(mesh *M), index sed_dupl (sed *A)
   
   	};
	
	/* GeMatrix */
    /*!  \class GeMatrix util_matrix.hpp "Include/util_matrix.hpp"
     *   \brief Matrix class for general dense matrixes
     *
     *	 The GeMatrix class is a general storage format for dense matrixes which stores
     *   all values including zero values. The storage format can be row or column major
     */	
    class GeMatrix {
    private:
    	long m; 			/*!< Number of Rows */
    	long n; 			/*!< Number of Columns */
    	long incRow;		/*!< Row Stride */
    	long incCol;		/*!< Column Stride */
        double *data; 		/*!< Numerical values */

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
    
    	/*!
    	 *   Constructor for GeMatrix
         *
         *   \param m Number of rows
         *   \param n Number of columns
         *   \param order Storage order ROWMAJOR | COLMAJOR
    	 */
    	GeMatrix(long m, long n, StorageOrder order) :
    			m(m), n(n), 
    			incRow(order == StorageOrder::COLMAJOR? 1: n),
    			incCol(order == StorageOrder::ROWMAJOR? 1: m),
    			data(new double[m*n]) {}
    	
    	/*!
    	 *   Constructor for GeMatrix from SedMatrix (always as colmajor)
         *
         *   \param sed Reference to SedMatrix
         *	 \param sym Is given SedMatrix symmetric
    	 */		
    	GeMatrix(SedMatrix &sed, bool sym) :
    			m(sed.get_n()), n(sed.get_n()),
    			incRow(1), incCol(n),
    			data(new double[sed.get_n() * sed.get_n()]) {
    		// Diagonal
    		for (long i = 0; i < n; ++i){
        		(*this)(i, i) = sed(i);
    		}
    		
    		// Super-/Sub-Diagonal
    		if (sym){
        		for (long j=0; j < n; ++j){
            		for (long p=sed.get_i(j); p < sed.get_i(j+1); ++p){
                		(*this)(sed.get_i(p), j) = sed(p);
                		(*this)(j, sed.get_i(p)) = sed(p);
            		}
        		}
    		} else {
        		for (long j=0; j < n; ++j){
            		for (long p=sed.get_i(j); p < sed.get_i(j+1); ++p){
                		(*this)(sed.get_i(p), j) = sed(p);
            		}
        		}
    		}			
    			
		}
    	
    	// Debug purpose
    	void Init() {
    		for (long i = 0; i < m; ++i) {
		     	for (long j = 0; j < n; ++j) {
		        	(*this)(i, j) = j + i * n;
		     	}
         	}
    	}	
        
        void Print();
   	};
}

#endif //HPC2_UTIL_MATRIX_HPP
