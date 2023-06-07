#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 double omega, long max_it, double tol) {

        long n = K.get_n();
        assert(n == f.count());

        // Declare
        BlasVector u(n);
        BlasVector r(n);

        // Calculate inverse diagonal
        BlasVector d = K.Diag();
        for (auto &value : d)
            value = 1 / value;

        // Initialise
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);

        for (long k = 1; k <= max_it; ++k) {
            if (r.Amax() < tol)
                break;

            // u = u + omega * d * r
            for (long i = 0; i < n; ++i)
                u(i) += omega * d(i) * r(i);
            // r = f - K * u
            r.Copy(f);
            K.SymSpmv(-1, u, 1, r);
        }

        return u;
    }
}