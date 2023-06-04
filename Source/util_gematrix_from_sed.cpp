#include "hpc.hpp"

namespace Util {
	
	/*
	 * Initialize Matrix from Sed Matrix
     *
     * sed: Reference to SedMatrix
     * sym: Is SedMatrix symmetric
     *
	 */
	void GeMatrix::from_sed(SedMatrix &sed, bool sym) {
		// Diagonal
		for (long i = 0; i < n; ++i){
    		(*this)(i, i) = sed(i);
		}
		
		// Super-/Sub-Diagonal
		if (sym){
    		for (long j=0; j < n; ++j){
        		for (long p=sed.get_ptr(j); p < sed.get_ptr(j+1); ++p){
            		(*this)(sed.get_ptr(p), j) = sed(p);
            		(*this)(j, sed.get_ptr(p)) = sed(p);
        		}
    		}
		} else {
    		for (long j=0; j < n; ++j){
        		for (long p=sed.get_ptr(j); p < sed.get_ptr(j+1); ++p){
            		(*this)(sed.get_ptr(p), j) = sed(p);
        		}
    		}
		}
	}  			
}
