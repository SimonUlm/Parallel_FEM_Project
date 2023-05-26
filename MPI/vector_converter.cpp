#include <cassert>

#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

namespace Conversion {

    void VectorConverter::AccumulatedToDistributed(std::vector<long> &local_vector) {
        assert(vector.size() == local_to_global_->count);

        for (long i = 0; i < local_vector.size(); ++i)
            local_vector[i] /= local_nodes_priority_(i);
    }

    void VectorConverter::DistributedToAccumulated(std::vector<long> &local_vector, MPI_Comm comm) {
        assert(vector.size() == local_to_global_->count);

        std::vector<long> global_vector_send(n_global_nodes_);
        for (long i = 0; i < local_vector.size(); ++i)
            global_vector_send[(*local_to_global_)(i)] = local_vector[i];

        std::vector<long> global_vector_recv(n_global_nodes_);
        MPI_Allreduce(global_vector_send.data(), global_vector_recv.data(), (int) n_global_nodes_,
                      MPI_LONG, MPI_SUM, comm);

        for (long i = 0; i < local_vector.size(); ++i)
            local_vector[i] = global_vector_send[(*local_to_global_)(i)];
    }
}
#endif // _MPI