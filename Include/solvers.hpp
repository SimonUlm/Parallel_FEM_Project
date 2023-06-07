#ifndef HPC2_SOLVERS_HPP
#define HPC2_SOLVERS_HPP

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>
#endif

namespace Solver {

    constexpr long kMaxItCG = 1e3;
    constexpr double kTolCG = 1e-5;
    constexpr double kOmegaJacobi = 1;
    constexpr long kMaxItJacobi = 1e6;
    constexpr double kTolJacobi = 1e-5;

    Util::BlasVector SolveCG(Util::SedMatrix &A, Util::BlasVector &r,
                             long max_it = kMaxItCG, double tol = kTolCG);
    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 double omega = kOmegaJacobi, long max_it = kMaxItJacobi, double tol = kTolJacobi);

#ifdef _MPI
    double ParallelDot(Util::Vector<double> &v_acc, Util::Vector<double> &v_dist);

    Util::BlasVector SolveCGParallel(Util::SedMatrix &K, Util::BlasVector &r,
                                     Skeleton::Skeleton &local_skel,
                                     long max_it = kMaxItCG, double tol = kTolCG);
    Util::BlasVector SolveJacobiParallel(Util::SedMatrix &K, Util::BlasVector &f,
                                         Skeleton::Skeleton &local_skel,
                                         double omega = kOmegaJacobi, long max_it = kMaxItJacobi, double tol = kTolJacobi);
#endif

}
#endif //HPC2_SOLVERS_HPP
