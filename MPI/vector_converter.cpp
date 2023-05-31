#include <cassert>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Mesh {

    void VectorConverter::AccumulatedToDistributed(Util::Vector<double> &local_vector) const {
        assert(local_vector.count() == local_to_global_->count());

        for (long i = 0; i < local_vector.count(); ++i)
            local_vector(i) /= (double) local_nodes_priority_(i);
    }

    void VectorConverter::DistributedToAccumulated(Util::Vector<double> &local_vector,
                                                   MPI_Comm comm, Skeleton::Skeleton &local_skel) const {
        assert(local_vector.count() == local_to_global_->count());

        // Create global crosspoint vector from local vector
        Util::Vector<double> global_vector_send(n_global_crosspoints_);
        for (auto point : local_skel.get_crosspoints())
            global_vector_send((*local_to_global_)(point)) = local_vector(point);

        // Allreduce
        Util::Vector<double> global_vector_recv(n_global_crosspoints_);
        MPI_Allreduce(global_vector_send.data(), global_vector_recv.data(),
                      (int) n_global_crosspoints_, MPI_DOUBLE, MPI_SUM, comm);

        // Write entries from global crosspoint vector back into local vector
        for (auto point : local_skel.get_crosspoints())
            local_vector(point) = global_vector_recv((*local_to_global_)(point));
    }
}
#endif // _MPI
