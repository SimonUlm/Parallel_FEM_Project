#include <algorithm>        // for std::copy
#include <cmath>            // for nan(), fabs()

#include "hpc.hpp"

namespace Util {
	
	// ----- Level 1 -----
	
    // x' * x
    double BlasVector::Dot(BlasVector &y) {
#ifndef NDEBUG
    	assert(count() == y.count());
#endif
    	
    	double sum = 0;
    	for (long i = 0; i < count(); ++i) {
    		sum += (*this)(i) * y(i);
    	}
    	
    	return sum;
    }
    
    // x <- alpha * x
    void BlasVector::Scal(double alpha) {
    	if (alpha == 0) {
    		for (long i = 0; i < count(); ++i) {
    			data_[i] = 0;
    		}
    	} else {
    		for (long i = 0; i < count(); ++i) {
    			data_[i] *= alpha;
    		}
    	}
    }
    
    // y <- alpha * x + y
    void BlasVector::Axpy(double alpha, BlasVector &x) {
    	if (alpha == 0){
    		return;
		}
		for (long i = 0; i < count(); ++i) {
			data_[i] += alpha*x(i);
		}
    
    }
    
    // max(x)
    double BlasVector::Amax() {
    	double max = 0;
    	
    	for (long i = 0; i < count(); ++i) {
			if (fabs(data_[i])>max) {
			    max = fabs(data_[i]);
			}
		}
		return max;	
    }
}
