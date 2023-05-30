#include <iostream>
#include <cstdio>

#include "hpc.hpp"

namespace Skeleton {

    void Skeleton::CreateLocal(int rank, Mesh::LocalMesh &local_mesh) {
        // Determine size of local skeleton
        long n_borders_loc = 0;

        for (long i = 0; i < n_borders; ++i) {
            if ((com_borders(i).get_L() == rank) || (com_borders(i).get_R() == rank)) {
                n_borders_loc += 1;
            }
        }

        // Allocate local skeleton
        long nodes_per_border = com_border_nodes.get_n_nodes();
        Skeleton local_skel(n_borders_loc, nodes_per_border, LOCAL);

        // Loop again to generate entries in local_skel
        long local_border_ix = 0;
        long global_node_ix = 0;
        long local_node_ix = 0;
        Util::List<long> &local2global = local_mesh.local_to_global;
        long length_l2g = local2global.count;

        for (long i = 0; i < n_borders; ++i) {
            if (com_borders(i).get_L() != rank && com_borders(i).get_R() != rank)
                continue;

            copy_border_entries(i, local_skel.get_border(local_border_ix));

            // Map c1 and c2 from global to local
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
                local_skel.set_border_node(local_border_ix, array_ix, local_node_ix);
            }
            local_border_ix += 1;
        }

        // Write local skeleton to current object
        // Global skeleton locally obsolete
        *this = std::move(local_skel);
    }
} // Namespace Skeleton

