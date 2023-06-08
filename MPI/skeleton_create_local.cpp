#include <cstdio>
#include <iostream>
#include <mpi.h>

#include "hpc.hpp"


namespace Skeleton {

    /*  
	 * Transforms global Skeleton to local Skeleton
     *
     * Creates a temporary Skeleton copies the ComBorder's which correspond with the local process
     * to the temporary Skeleton. At the end the old global Skeleton is overwritten with the local one.
     *
     * local_mesh: LocalMesh to map global and local indexing
     * 
	 */
    void Skeleton::CreateLocal(Mesh::LocalMesh &local_mesh) {
        // Determine size of local skeleton (and remember cross points)
        long n_borders_loc = 0;
        for (long i = 0; i < n_borders_; ++i) {
            if ((com_borders_(i).L() == rank_) || (com_borders_(i).R() == rank_)) {
                n_borders_loc += 1;
                crosspoints_(com_borders_(i).c1()) = true;
                crosspoints_(com_borders_(i).c2()) = true;
            }
        }

        // Allocate local skeleton
        long nodes_per_border = com_border_nodes_.n_nodes();

        Skeleton local_skel(n_borders_loc, nodes_per_border, comm_, rank_, LOCAL);

        // Loop again to generate entries in local_skel
        long local_border_ix = 0;
        long global_node_ix = 0;
        long local_node_ix = 0;
        Util::Vector<long> &local2global = local_mesh.local_to_global;
        long length_l2g = local2global.count();

        for (long i = 0; i < n_borders_; ++i) {
            if (com_borders_(i).L() != rank_ && com_borders_(i).R() != rank_)
                continue;

            // Corresponding ComBorder Found -> Copy entries to new local skeleton
            com_borders_(i).copy_entries(local_skel.com_borders_(local_border_ix));

            // Map c1_ and c2_ from global to local and set index_ for local skeleton
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == com_borders_(i).c1()) {
                    local_skel.com_borders_(local_border_ix).set_c1(j);
                    break;
                }
            }
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == com_borders_(i).c2()) {
                    local_skel.com_borders_(local_border_ix).set_c2(j);
                    break;
                }
            }
            local_skel.com_borders_(local_border_ix).set_index(local_border_ix);

            // Copy border nodes_ and map from global to local node index_
            for (long array_ix = 0; array_ix < nodes_per_border; ++array_ix) {
                global_node_ix = com_border_nodes_.get_entry(i, array_ix);
                for (long j = 0; j < length_l2g; ++j) {
                    if (local2global(j) == global_node_ix) {
                        local_node_ix = j;
                        break;
                    }
                    if (j == length_l2g - 1) {
                        std::cerr << "Skeleton.CreateLocal(...): MAPPING FAILED\n" << std::endl;
                        abort();
                    }
                }
                // Write local node ix to local skel
                local_skel.com_border_nodes_.set_entry(local_border_ix, array_ix, local_node_ix);
            }
            local_border_ix += 1;
        }

        // Count cross points
        long count = 0;
        for (auto &point: crosspoints_)
            if (point == true)
                ++count;

        // Create local cross points
        local_skel.crosspoints_ = Util::Vector<long>(count);
        long crosspoint_idx = 0;
        for (long i = 0; i < crosspoints_.count(); ++i) {
            if (crosspoints_(i) == false)
                continue;
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == i) {
                    local_skel.crosspoints_(crosspoint_idx) = j;
                    ++crosspoint_idx;
                    break;
                }
            }
        }

        // Write local skeleton to current object
        // Global skeleton locally obsolete
        *this = std::move(local_skel);
    }
}


