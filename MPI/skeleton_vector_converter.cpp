#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

namespace Skeleton {

    void VectorConverter::AccumulatedToDistributed(Util::Vector<double> &local_vector) const {
#ifndef NDEBUG
        assert(local_vector.count() == local_to_global_->count());
#endif

        for (long i = 0; i < local_vector.count(); ++i)
            local_vector(i) /= (double) local_nodes_priority_(i);
    }

    void VectorConverter::DistributedToAccumulated(Util::Vector<double> &local_vector,
                                                   const Skeleton &local_skel) const {
        DistributedToAccumulated(local_vector, local_vector, local_skel);
    }

    void VectorConverter::DistributedToAccumulated(Util::Vector<double> &local_vector_send,
                                                   Util::Vector<double> &local_vector_recv,
                                                   const Skeleton &local_skel) const {
#ifndef NDEBUG
        assert(local_vector_send.count() == local_to_global_->count());
        assert(local_vector_recv.count() == local_to_global_->count());
#endif

        MPI_Comm comm = local_skel.comm();
        int rank = local_skel.rank();

        // Copy all regular nodes_ from send vector to recv vector
        local_vector_recv.Copy(local_vector_send);

        // Create global cross point vector from local vector
        Util::Vector<double> global_vector_send(n_global_crosspoints_);
        for (auto &point: local_skel.crosspoints())
            global_vector_send((*local_to_global_)(point)) = local_vector_send(point);

        // Allreduce
        Util::Vector<double> global_vector_recv(n_global_crosspoints_);
        MPI_Allreduce(global_vector_send.data(), global_vector_recv.data(),
                      (int) n_global_crosspoints_, MPI_DOUBLE, MPI_SUM, comm);

        // Write entries from global cross point vector back into local vector
        for (auto &point: local_skel.crosspoints())
            local_vector_recv(point) = global_vector_recv((*local_to_global_)(point));

        // Exchange interface points with neighbors
        MPI_Status status;
        long n_border_nodes = local_skel.n_border_nodes();
        Util::Vector<double> border_nodes_vector_send(n_border_nodes);
        Util::Vector<double> border_nodes_vector_recv(n_border_nodes);

        for (long col = 0; col < kNumberOfColors; ++col) {
            // Find border_ix and neighbor to exchange with
            // Initialise neighbor with own rank
            int neighbor = rank;
            long border_ix;
            for (long i = 0; i < local_skel.n_borders(); ++i) {
                auto &border = local_skel.get_border(i);
                if (border.color() != col)
                    continue;
                neighbor = (int) ((border.L() == rank) ? border.R() : border.L());
                border_ix = i;
                break;
            }

            // Check whether border exits
            if (neighbor == rank)
                continue;

            // Exchange data
            for (long i = 0; i < n_border_nodes; ++i)
                border_nodes_vector_send(i) = local_vector_send(local_skel.get_border_node(border_ix, i));
            MPI_Sendrecv(border_nodes_vector_send.data(), (int) n_border_nodes, MPI_DOUBLE, neighbor, 0,
                         border_nodes_vector_recv.data(), (int) n_border_nodes, MPI_DOUBLE, neighbor, 0, comm, &status);
            for (long i = 0; i < n_border_nodes; ++i)
                local_vector_recv(local_skel.get_border_node(border_ix, i)) += border_nodes_vector_recv(i);
        }
    }

    void VectorConverter::GatherAccumulatedVector(Util::Vector<double> &local_vector_send,
                                                  Util::Vector<double> &global_vector_recv,
                                                  const Skeleton &local_skel) const {
#ifndef NDEBUG
        assert(local_vector_send.count() == local_to_global_->count());
#endif

        // Convert vector to ease the process
        AccumulatedToDistributed(local_vector_send);

        GatherDistributedVector(local_vector_send, global_vector_recv, local_skel);

    }

    void VectorConverter::GatherDistributedVector(Util::Vector<double> &local_vector_send,
                                                  Util::Vector<double> &global_vector_recv,
                                                  const Skeleton &local_skel) const {
#ifndef NDEBUG
        assert(local_vector_send.count() == local_to_global_->count());
#endif

        MPI_Comm comm = local_skel.comm();
        int rank = local_skel.rank();

        // Create global vector from local vector
        Util::Vector<double> global_vector_send(n_global_nodes_);
        for (long i = 0; i < local_vector_send.count(); ++i)
            global_vector_send((*local_to_global_)(i)) = local_vector_send(i);

        // Allreduce
        if (rank == 0)
            global_vector_recv = Util::Vector<double>(n_global_nodes_);
        MPI_Reduce(global_vector_send.data(), global_vector_recv.data(), (int) n_global_nodes_,
                   MPI_DOUBLE, MPI_SUM, 0, comm);
    }
}
#endif // _MPI
