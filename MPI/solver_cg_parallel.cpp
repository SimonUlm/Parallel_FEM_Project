#include <cassert>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Util;

namespace Solver {

    BlasVector SolveCGParallel(SedMatrix &K, BlasVector &f,
                               MPI_Comm comm, int rank,
                               Mesh::LocalMesh &local_mesh, Skeleton::Skeleton &local_skel,
                               long max_it, double tol) {

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
        BlasVector w(n);
        BlasVector s(n);
        BlasVector v(n);
        double sigma;
        double old_sigma;
        double alpha;

        // Initialise
        // r = f - K * u
        r.Copy(f);
        K.SymSpmv(-1, u, 1, r);
        // w = r
        local_mesh.vector_converter().DistributedToAccumulated(r, w, comm, rank, local_skel);
        // s = w
        s.Copy(w);
        // sigma = <w, r>
        sigma = w.Dot(r);

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
            // w = r
            local_mesh.vector_converter().DistributedToAccumulated(r, w, comm, rank, local_skel);
            // sigma = <w, r>
            old_sigma = sigma;
            sigma = w.Dot(r);
            // s = r + (sigma / old_sigma) * s
            s.Scal(sigma / old_sigma);
            s.Axpy(1, w);
        }

        return u;
    }
}

#endif