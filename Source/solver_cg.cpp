#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    BlasVector SolveCG(SedMatrix &K, BlasVector &f, long max_it, double tol) {

        constexpr double kTol = 1e-5;

        // Handle input
        long n = K.get_n();
        assert(n == f.count());
        if (max_it == 0)
            max_it = n;
        if (tol == 0)
            tol = kTol;

        // Declare
        BlasVector u(n);
        BlasVector r(n);
        BlasVector s(n);
        BlasVector v(n);
        double sigma;
        double old_sigma;
        double alpha;

        // Initialise
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);
        // s = r
        s.Copy(r);
        // sigma = <r, r>
        sigma = r.Dot(r);

        for (long k = 1; k <= max_it; ++k) {
            if (r.Amax() < tol)
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
            old_sigma = sigma;
            sigma = r.Dot(r);
            // s = r + (sigma / old_sigma) * s
            s.Scal(sigma / old_sigma);
            s.Axpy(1, r);
        }

        return std::move(u);
    }
}