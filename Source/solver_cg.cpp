#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    BlasVector SolveCG(SedMatrix &A, BlasVector &r, long max_it, double tol) {

        long n = A.get_n();
        assert(n == r.count());

        // Declare
        BlasVector x(n);
        BlasVector q(n);
        BlasVector p(n);
        BlasVector v(n);
        double sigma;
        double sigma_old;
        double sigma_0;
        double alpha;

        // Initialise
        // q = r - A * x
        q.Copy(r);
        A.SymSpmv(-1, x, 1, q);
        // p = q
        p.Copy(q);
        // sigma = <q, q>
        sigma = q.Dot(q);
        sigma_old = sigma;
        sigma_0 = sigma;

        for (long k = 1; k <= max_it; ++k) {
            if (std::sqrt(sigma / sigma_0) < tol)
                break;

            // v = A * p
            A.SymSpmv(1, p, 0, v);
            // alpha = sigma / <p, v>
            alpha = sigma / p.Dot(v);
            // x = x + alpha * p;
            x.Axpy(alpha, p);

            // q = q - alpha * v
            q.Axpy(-alpha, v);
            // sigma = <q, q>
            sigma = q.Dot(q);
            // p = q + (sigma / sigma_old) * p
            p.Scal(sigma / sigma_old);
            p.Axpy(1, q);
            // sigma_old = sigma
            sigma_old = sigma;
        }

        return x;
    }
}