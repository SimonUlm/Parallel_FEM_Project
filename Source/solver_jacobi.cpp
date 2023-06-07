#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 double &error,
                                 double omega, long max_it, double tol) {

        long n = K.get_n();
#ifndef NDEBUG
        assert(n == f.count());
#endif

        // Declare
        BlasVector r(n); // residuum
        BlasVector u(n); // solution
        BlasVector d(n); // diagonal of the stiffness matrix
        double sigma;

        // Initialise
        // d = diag(K)^-1
        d = K.Diag();
        for (auto &value : d)
            value = 1 / value;
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);
        // sigma = <r, r>
        sigma = r.Dot(r);

        // Iterate
        for (long k = 1; k <= max_it; ++k) {
            if (sigma < tol * tol)
                break;

            // u = u + omega * d * r
            for (long i = 0; i < n; ++i)
                u(i) += omega * d(i) * r(i);

            // r = f - K * u
            r.Copy(f);
            K.SymSpmv(-1, u, 1, r);
            // sigma = <r, r>
            sigma = r.Dot(r);
        }

        error = sigma;
        return u;
    }
}