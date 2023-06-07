#include <cassert>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Util;

namespace Solver {

    BlasVector SolveCgParallel(Util::SedMatrix &K, Util::BlasVector &r,
                               Skeleton::Skeleton &local_skel,
                               double &error,
                               long max_it, double tol) {

        long n = K.get_n();
#ifndef NDEBUG
        assert(n == r.count());
#endif

        // Declare
        // The input vector r (namely the rhs) will be used as residuum
        BlasVector w(n); // accumulated residuum
        BlasVector u(n); // solution
        BlasVector s(n); // direction
        BlasVector v(n); // temporary vector for storing K * u
        double sigma;
        double sigma_old;
        double alpha;    // step size

        // Initialise
        // r = r - K * u
        K.SymSpmv(-1, u, 1, r);
        // w = r
        w.Copy(r);
        local_skel.DistributedToAccumulated(w);
        // s = w
        s.Copy(w);
        // sigma = <w, w>
        sigma = Solver::ParallelDot(w, r);

        // Iterate
        for (long k = 1; k <= max_it; ++k) {
            if (sigma < tol * tol)
                break;

            // v = K * s
            K.SymSpmv(1, s, 0, v);
            // alpha = sigma / <s, v>
            alpha = sigma / Solver::ParallelDot(s, v);
            // u = u + alpha * s;
            u.Axpy(alpha, s);

            // r = r - alpha * v
            r.Axpy(-alpha, v);
            // w = r
            w.Copy(r);
            local_skel.DistributedToAccumulated(w);
            // sigma = <w, r>
            sigma_old = sigma;
            sigma = Solver::ParallelDot(w, r);
            // s = w + (sigma / sigma_old) * s
            s.Scal(sigma / sigma_old);
            s.Axpy(1, w);
        }

        error = sigma;
        return u;
    }
}

#endif