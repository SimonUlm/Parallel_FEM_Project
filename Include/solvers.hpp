#ifndef HPC2_SOLVERS_HPP
#define HPC2_SOLVERS_HPP

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>
#endif

namespace Solver {

    Util::BlasVector SolveCG(Util::SedMatrix &K, Util::BlasVector &f,
                             long max_it = 0, double tol = 0);
    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 double omega = 0, long max_it = 0, double tol = 0);

#ifdef _MPI
    Util::BlasVector SolveCGParallel(Util::SedMatrix &K, Util::BlasVector &f,
                                     Skeleton::Skeleton &local_skel,
                                     long max_it, double tol);
    Util::BlasVector SolveJacobi(Util::SedMatrix &K, Util::BlasVector &f,
                                 Skeleton::Skeleton &local_skel,
                                 double omega = 0, long max_it = 0, double tol = 0);
#endif

}
#endif //HPC2_SOLVERS_HPP
