#include <cstdio>

#include "hpc.hpp"

namespace Util {	
	/*
     *   Prints Data of GeMatrix
     */
    void BlasVector::Print() {
        // Start printing
        printf("\n=========== Print BlasVector Data ===========\n");
        
        // Print Matrix sizes
        printf("Number of entries: %zu\n", count_);
               
        // Print Data
        printf("\nData:\n");
        for (long i = 0; i < count_; ++i) {
		 	printf("  %.4d \n", (*this)(i));
         }
        
        /* 
        // Print overall storage requirements
        printf("\nMemory\n");
        printf("Additional info	: %12zu Byte\n", 4 * sizeof(long));
        printf("Data		: %12zu Byte\n", m * n * sizeof(double));
        long total = 4 * sizeof(long) + m * n * sizeof(double);
        printf("Total       	: %12.6g MByte\n", (double) total / 1024. / 1024.);
        */
    }
}
