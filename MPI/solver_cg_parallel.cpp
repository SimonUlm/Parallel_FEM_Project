#include <cassert>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Util;

namespace Solver {

    BlasVector SolveCGParallel(SedMatrix &K, BlasVector &r,
                               Skeleton::Skeleton &local_skel,
                               long max_it, double tol) {

        long n = K.get_n();
        assert(n == r.count());

        // Declare
        BlasVector x(n);
        BlasVector q(n);
        //BlasVector w(n);
        BlasVector p(n);
        BlasVector a(n);
        double rho_new;
        double rho_old;
        double lambda;

        // Initialise
        // q = r - K * x
        K.SymSpmv(1, x, 0, a);
        r.Axpy(-1, a);
        q.Copy(r);
        // w = q
        local_skel.DistributedToAccumulated(q);
        // p = w
        p.Copy(q);
        // rho_new = <w, q>
        rho_new = Solver::ParallelDot(q, r);

        for (long k = 1; k <= max_it; ++k) {
            if (rho_new < tol * tol)
                break;

            // a = K * p
            K.SymSpmv(1, p, 0, a);
            // lambda = rho_new / <p, a>
            lambda = rho_new / Solver::ParallelDot(p, a);
            // x = x + lambda * p;
            x.Axpy(lambda, p);
            // q = q - lambda * a
            r.Axpy(-lambda, a);
            // w = q
            q.Copy(r);
            local_skel.DistributedToAccumulated(q);
            // rho_new = <w, q>
            rho_old = rho_new;
            rho_new = Solver::ParallelDot(q, r);
            // p = q + (rho_new / rho_old) * p
            p.Scal(rho_new / rho_old);
            p.Axpy(1, q);
        }

        return x;
    }
}

#endif