#ifndef HPC2_SOLVERS_HPP
#define HPC2_SOLVERS_HPP

#include "hpc.hpp"

namespace Solver {

    Util::BlasVector SolveCG(Util::SedMatrix &K, Util::BlasVector &f, long max_it = 0, double tol = 0);
}

#endif //HPC2_SOLVERS_HPP
