#include "hpc.hpp"

#define N 1000
#define TOL 0.00001

namespace Util {
	
	/*
	 * Solve SedMatrix with CG
     *
     * rhs: Right hand side of equation
     * sol: Pointer to solution vector with start value when calling function
     *
	 */
	void SedMatrix::SolveCg(BlasVector &rhs, BlasVector &sol) {
		
		BlasVector r(rhs.length());
		BlasVector d(rhs.length());
		BlasVector a(rhs.length());
		
		// r <- rhs
		r.Copy(rhs);
		
		// d <- r
		d.Copy(r);
		
		// r <- - this * x_0 + r
		SymSpmv(-1, sol, 1, r);	
		
		double rho = r.Dot(r);
		double alpha = 0;
		double rho_new = 0;
		
		for (long i = 0; i < N; ++i) {
			if (r.Amax() < TOL) {
				break;
			}
			
			// a <- this * d
			SymSpmv(1, d, 0, a);
			
			alpha = rho / d.Dot(a);
			
			// sol <- alpha * d + sol
			sol.Axpy(alpha, d);
			
			// r <- -alpha * a + r
			r.Axpy(-1*alpha, a);
			
			rho_new = r.Dot(r);
			
			// d <- r + rho_new/rho * d
			d.Scal(rho_new/rho);
			d.Axpy(1, r);
			
			rho = rho_new;
		}
	}
}
