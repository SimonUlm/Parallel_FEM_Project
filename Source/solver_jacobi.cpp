#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 double omega, long max_it, double tol) {

        constexpr double kOmega = 0.2;
        constexpr long kMaxIt = 1000;
        constexpr double kTol = 1e-5;

        // Handle input
        long n = K.get_n();
        assert(n == f.count());
        if (omega == 0)
            omega = kOmega;
        if (max_it == 0)
            max_it = kMaxIt;
        if (tol == 0)
            tol = kTol;

        // Declare
        BlasVector u(n);
        BlasVector r(n);

        // Calculate inverse diagonal
        BlasVector d = K.Diag();
        for (auto &value : d)
            {value = 1 / value; printf("%lf\n", value);}

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