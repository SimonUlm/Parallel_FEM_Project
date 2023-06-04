#ifndef HPC2_SOLVERS_HPP
#define HPC2_SOLVERS_HPP

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>
#endif

namespace Solver {

    Util::BlasVector SolveCG(Util::SedMatrix &K, Util::BlasVector &f, long max_it = 0, double tol = 0);

#ifdef _MPI
    Util::BlasVector SolveCGParallel(Util::SedMatrix &K, Util::BlasVector &f,
                                     Mesh::LocalMesh &local_mesh, Skeleton::Skeleton &local_skel,
                                     long max_it, double tol);
#endif

}
#endif //HPC2_SOLVERS_HPP
