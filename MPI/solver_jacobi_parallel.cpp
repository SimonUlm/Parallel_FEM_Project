#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

using namespace Util;

namespace Solver {

    Util::BlasVector SolveJacobiParallel(Util::SedMatrix &K, Util::BlasVector &f,
                                         Skeleton::Skeleton &local_skel,
                                         double &error,
                                         double omega, long max_it, double tol) {

        long n = K.n();
#ifndef NDEBUG
        assert(n == f.count());
#endif

        // Declare
        BlasVector r(n); // residuum
        BlasVector w(n); // accumulated residuum
        BlasVector u(n); // solution
        BlasVector d(n); // diagonal of the stiffness matrix
        double sigma;

        // Initialise
        // d = diag(K)^-1
        d = K.Diag();
        local_skel.DistributedToAccumulated(d);
        for (auto &value: d)
            value = 1 / value;
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);
        // w = r
        local_skel.DistributedToAccumulated(r, w);
        // sigma = <w, r>
        sigma = Solver::ParallelDot(w, r, local_skel);

        // Iterate
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
            sigma = Solver::ParallelDot(w, r, local_skel);
        }

        error = sigma;
        return u;
    }
}

#endif