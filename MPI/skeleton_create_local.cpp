#include <cstdio>
#include <iostream>
#include <mpi.h>

#include "hpc.hpp"


namespace Skeleton {

    /*  
	 * Transforms global Skeleton to local Skeleton
     *
     * Creates a tempory Skeleton copies the ComBorder's which correspond with the local process
     * to the tempory Skeleton. At the end the old global Skeleton is overwritten with the local one.
     *
     * rank: Number of Local Process
     * local_mesh: LocalMesh to map global and local indexing
     * 
	 */
    void Skeleton::CreateLocal(Mesh::LocalMesh &local_mesh) {
        // Determine size of local skeleton (and remember cross points)
        long n_borders_loc = 0;
        for (long i = 0; i < n_borders; ++i) {
            if ((com_borders(i).get_L() == rank) || (com_borders(i).get_R() == rank)) {
                n_borders_loc += 1;
                crosspoints(com_borders(i).get_c1()) = true;
                crosspoints(com_borders(i).get_c2()) = true;
            }
        }

        // Allocate local skeleton
        long nodes_per_border = com_border_nodes.get_n_nodes();

        Skeleton local_skel(n_borders_loc, nodes_per_border, comm, rank, LOCAL);

        // Loop again to generate entries in local_skel
        long local_border_ix = 0;
        long global_node_ix = 0;
        long local_node_ix = 0;
        Util::Vector<long> &local2global = local_mesh.local_to_global;
        long length_l2g = local2global.count();

        for (long i = 0; i < n_borders; ++i) {
            if (com_borders(i).get_L() != rank && com_borders(i).get_R() != rank)
                continue;
			
			// Corresponding ComBorder Found -> Copy entries to new local skeleton
            com_borders(i).copy_entries(local_skel.com_borders(local_border_ix));

            // Map c1 and c2 from global to local and set index for local skeleton
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == com_borders(i).get_c1()) {
                    local_skel.com_borders(local_border_ix).set_c1(j);
                    break;
                }
            }
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == com_borders(i).get_c2()) {
                    local_skel.com_borders(local_border_ix).set_c2(j);
                    break;
                }
            }
            // Refactor: index variable actually unnecessary
            local_skel.com_borders(local_border_ix).set_index(local_border_ix);

            // Copy border nodes and map from global to local node index
            for (long array_ix = 0; array_ix < nodes_per_border; ++array_ix) {
                global_node_ix = com_border_nodes.get_entry(i, array_ix);
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
                local_skel.com_border_nodes.set_entry(local_border_ix, array_ix, local_node_ix);
            }
            local_border_ix += 1;
        }

        // Count crosspoints
        long count = 0;
        for (auto &point : crosspoints)
            if (point == true)
                ++count;

        // Create local crosspoints
        local_skel.crosspoints = Util::Vector<long>(count);
        long crosspoint_idx = 0;
        for (long i = 0; i < crosspoints.count(); ++i) {
            if (crosspoints(i) == false)
                continue;
            for (long j = 0; j < length_l2g; ++j) {
                if (local2global(j) == i) {
                    local_skel.crosspoints(crosspoint_idx) = j;
                    ++crosspoint_idx;
                    break;
                }
            }
        }

        // Write local skeleton to current object
        // Global skeleton locally obsolete
        *this = std::move(local_skel);
    }
} // Namespace Skeleton


