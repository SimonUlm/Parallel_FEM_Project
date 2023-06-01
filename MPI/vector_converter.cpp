#include <cassert>
#include <cstdio>

#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

namespace Mesh {

    void VectorConverter::AccumulatedToDistributed(Util::Vector<double> &local_vector) const {
        assert(local_vector.count() == local_to_global_->count());

        for (long i = 0; i < local_vector.count(); ++i)
            local_vector(i) /= (double) local_nodes_priority_(i);
    }

    void VectorConverter::DistributedToAccumulated(Util::Vector<double> &local_vector, MPI_Comm comm,
                                                   int rank, Skeleton::Skeleton &local_skel) const {
        assert(local_vector.count() == local_to_global_->count());

        // Create global cross point vector from local vector
        Util::Vector<double> global_vector_send(n_global_crosspoints_);
        for (auto &point : local_skel.get_crosspoints())
            global_vector_send((*local_to_global_)(point)) = local_vector(point);

        // Allreduce
        Util::Vector<double> global_vector_recv(n_global_crosspoints_);
        MPI_Allreduce(global_vector_send.data(), global_vector_recv.data(),
                      (int) n_global_crosspoints_, MPI_DOUBLE, MPI_SUM, comm);

        // Write entries from global cross point vector back into local vector
        for (auto &point : local_skel.get_crosspoints())
            local_vector(point) = global_vector_recv((*local_to_global_)(point));

        // Exchange interface points with neighbors
        MPI_Status status;
        long n_border_nodes = local_skel.get_n_border_nodes();
        Util::Vector<double> send_vector(n_border_nodes);
        Util::Vector<double> recv_vector(n_border_nodes);

        for (long col = 0; col < Skeleton::kNumberOfColors; ++col) {
            // Find border_ix and neighbor to exchange with
            // Initialise neighbor with own rank
            int neighbor = rank;
            long border_ix;
            for (long i = 0; i < local_skel.get_n_borders(); ++i) {
                auto &border = local_skel.get_border(i);
                if (border.get_color() != col)
                    continue;
                neighbor = (int) (border.get_L() == rank) ? border.get_R() : border.get_L();
                border_ix = i;
                break;
            }

            // Check whether border exits
            if (neighbor == rank)
                continue;

            // Exchange data
            for (long i = 0; i < n_border_nodes; ++i)
                send_vector(i) = local_vector(local_skel.get_border_node(border_ix, i));
            MPI_Sendrecv(send_vector.data(), (int) n_border_nodes, MPI_DOUBLE, neighbor, 0,
                         recv_vector.data(), (int) n_border_nodes, MPI_DOUBLE, neighbor, 0, comm, &status);
            for (long i = 0; i < n_border_nodes; ++i)
                local_vector(local_skel.get_border_node(border_ix, i)) += recv_vector(i);
        }
    }
}
#endif // _MPI
