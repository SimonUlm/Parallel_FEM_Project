#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    BlasVector SolveCg(Util::SedMatrix &K, Util::BlasVector &r, long max_it, double tol) {

        long n = K.get_n();
        assert(n == r.count());

        // Declare
        // The input vector r (namely the rhs) will be used as residuum
        BlasVector u(n); // solution
        BlasVector s(n); // direction
        BlasVector v(n); // temporary vector for storing K * u
        double sigma;
        double sigma_old;
        double alpha;    // step size

        // Initialise
        // r = r - K * u
        K.SymSpmv(-1, u, 1, r);
        // s = r
        s.Copy(r);
        // sigma = <r, r>
        sigma = r.Dot(r);

        // Iterate
        for (long k = 1; k <= max_it; ++k) {
            if (sigma < tol * tol)
                break;

            // v = K * s
            K.SymSpmv(1, s, 0, v);
            // alpha = sigma / <s, v>
            alpha = sigma / s.Dot(v);
            // u = u + alpha * s;
            u.Axpy(alpha, s);

            // r = r - alpha * v
            r.Axpy(-alpha, v);
            // sigma = <r, r>
            sigma_old = sigma;
            sigma = r.Dot(r);
            // s = r + (sigma / sigma_old) * s
            s.Scal(sigma / sigma_old);
            s.Axpy(1, r);
        }

        return u;
    }
}