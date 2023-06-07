#include <cassert>

#include "hpc.hpp"

using namespace Util;

namespace Solver {

    Util::BlasVector SolveJacobiParallel(Util::SedMatrix &K, Util::BlasVector &f,
                                         Skeleton::Skeleton &local_skel,
                                         double omega, long max_it, double tol) {

        long n = K.get_n();
        assert(n == f.count());

        // Declare
        BlasVector u(n);
        BlasVector r(n);
        BlasVector w(n);
        double sigma;

        // Calculate inverse diagonal and accumulate
        BlasVector d = K.Diag();
        local_skel.DistributedToAccumulated(d);
        for (auto &value : d)
            value = 1 / value;

        // Initialise
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);
        // w = r
        local_skel.DistributedToAccumulated(r, w);
        // sigma = <w, r>
        sigma = Solver::ParallelDot(w, r);

        for (long k = 1; k <= max_it; ++k) {
            if (sigma < tol * tol)
                break;

            // u = u + omega * d * r
            for (long i = 0; i < n; ++i)
                u(i) += omega * d(i) * w(i);
            // r = f - K * u
            r.Copy(f);
            K.SymSpmv(-1, u, 1, r);
            // w = r
            local_skel.DistributedToAccumulated(r, w);
            // sigma = <w, r>
            sigma = Solver::ParallelDot(w, r);
        }

        return u;
    }
}