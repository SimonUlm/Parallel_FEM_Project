#include <cassert>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Util {

    void VectorConverter::AccumulatedToDistributed(std::vector<double> &local_vector) const {
        assert(local_vector.size() == local_to_global_->count());

        for (long i = 0; i < local_vector.size(); ++i)
            local_vector[i] /= (double) local_nodes_priority_(i);
    }

    void VectorConverter::DistributedToAccumulated(std::vector<double> &local_vector, MPI_Comm comm) const {
        assert(local_vector.size() == local_to_global_->count());

        // Create global vector from local vector
        std::vector<double> global_vector_send(n_global_nodes_);
        for (long i = 0; i < local_vector.size(); ++i)
            global_vector_send[(*local_to_global_)(i)] = local_vector[i];

        // Allreduce
        std::vector<double> global_vector_recv(n_global_nodes_);
        MPI_Allreduce(global_vector_send.data(), global_vector_recv.data(), (int) n_global_nodes_,
                      MPI_DOUBLE, MPI_SUM, comm);

        // Write entries from global vector back into local vector
        for (long i = 0; i < local_vector.size(); ++i)
            local_vector[i] = global_vector_recv[(*local_to_global_)(i)];
    }
}
#endif // _MPI
