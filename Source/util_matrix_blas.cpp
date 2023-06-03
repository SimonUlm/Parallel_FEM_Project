#include <math.h>           // for nan(), fabs()

#include "hpc.hpp"

namespace Util {
	
	// ----- Level 1 -----
	
    // x' * x
    double BlasVector::Dot(BlasVector &y) {
    	assert(length() == y.length());
    	
    	double sum = 0;
    	for (long i = 0; i < length(); ++i) {
    		sum += (*this)(i) * y(i);
    	}
    	
    	return sum;
    }
    
    // x <- y
    void BlasVector::Copy(BlasVector &y) {
    	assert(length() == y.length());
    	for (long i = 0; i < y.length(); ++i) {
			data[i] = y(i);
		}
    }
    
    // x <- alpha * x
    void BlasVector::Scal(double alpha) {
    	if (alpha == 0) {
    		for (long i = 0; i < length(); ++i) {
    			data[i] = 0;
    		}
    	} else {
    		for (long i = 0; i < length(); ++i) {
    			data[i] *= alpha;
    		}
    	}
    }
    
    // y <- alpha * x + y
    void BlasVector::Axpy(double alpha, BlasVector &x) {
    	if (alpha == 0){
    		return;
		}
		for (long i = 0; i < length(); ++i) {
			data[i] += alpha*x(i);
		}
    
    }
    
    // max(x)
    double BlasVector::Amax() {
    	double max = 0;
    	
    	for (long i = 0; i < length(); ++i) {
			if (fabs(data[i])>max) {
			    max = fabs(data[i]);
			}
		}
		return max;	
    }
    
    
    // ----- Level 2 -----
    
    // y <- alpha * this * x + beta * y
    void SedMatrix::SymSpmv(double alpha, BlasVector &x, double beta, BlasVector &y) {
        
		y.Scal(beta);
		
		if (alpha == 0) {return;}
		
		for (long j = 0; j < n; ++j) {
			// Diagonal
			y(j) += alpha * data[j] * x(j);
			
			// Sub/Super-Diagonal
			for (long p = ptr_ind[j]; p < ptr_ind[j+1]; ++p) {
				// lower part axpy-based
				y(ptr_ind[p]) += alpha * data[p] * x(j);
				
				// upper part ddot-based
				y(j) += alpha * data[p] * x(ptr_ind[p]);
			}        		
		}
	}		

}
