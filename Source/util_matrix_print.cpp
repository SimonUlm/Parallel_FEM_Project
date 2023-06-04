#include <cstdio>

#include "hpc.hpp"

namespace Util {	
	/*
     *   Prints Data of GeMatrix
     */
    void GeMatrix::Print() {
        // Start printing
        printf("\n=========== Print GeMatrix Data ===========\n");
        
        // Print Matrix Storage Order
        printf("Storage Order:\t");
        printf(incRow == 1? "Column Major\n": "Row Major\n");
        
        // Print Matrix sizes
        printf("Matrix size:\n");
        printf("Number of rows: %zu\nNumber of columns: %zu\n", m, n);
               
        // Print Data
        printf("\nData:\n");
        for (long i = 0; i < m; ++i) {
         	printf("  ");
         	for (long j = 0; j < n; ++j) {
            	printf(" %4.5lf", (*this)(i, j));
         	}
         	printf("\n");
         }
         
        // Print overall storage requirements
        printf("\nMemory\n");
        printf("Additional info	: %12zu Byte\n", 4 * sizeof(long));
        printf("Data		: %12zu Byte\n", m * n * sizeof(double));
        long total = 4 * sizeof(long) + m * n * sizeof(double);
        printf("Total       	: %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
    
    /*
     *   Prints Data of SedMatrix
     */
    void SedMatrix::Print() {
    	// Start printing
        printf("\n=========== Print SedMatrix Data ===========\n");
        
        printf("Basic Information:\t %zu-by-%zu, nzmax: %zu nnz: %zu\n", n, n, nzmax, ptr_ind[n]);
        
        // Print diagonal
    	printf ("Diagonal entries:\n"); 
    	for (long j = 0 ; j < n ; j++) {
        	printf ("      %zu : %g\n", j, data ? data[j] : 1);
        }
        
        // Print off-diagonal
    	printf ("Off-diagonal entries:\n"); 
    	for (long j = 0 ; j < n ; j++) {
        	printf("    col %zu : locations %zu to %zu\n", j, ptr_ind[j], (ptr_ind[j+1]-1)) ;
        	for (long p = ptr_ind[j] ; p < ptr_ind[j+1] ; p++) {
            	printf("(%3zu,%3zu) : %g\n", ptr_ind[p], j, data[p]) ;
        	}
    	}
    	
    	// Print overall storage requirements
        printf("\nMemory\n");
        printf("Additional info			: %12zu Byte\n", 2 * sizeof(long));
        printf("Col Pointers and Row indices	: %12zu Byte\n", nzmax * sizeof(long));
        printf("Data				: %12zu Byte\n", nzmax * sizeof(double));
        long total = 2 * sizeof(long) + nzmax * sizeof(long) + nzmax * sizeof(double);
        printf("Total       			: %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
}
