#include "hpc.hpp"
#include <iostream>
#include <cstdio>

namespace Skeleton {
    void Skeleton::CreateLocal(long process, long* local2global, 
			       long length_l2g) {
	// Determine size of local skeleton
	long n_borders_loc = 0;

	for (long i = 0; i < n_borders; ++i) {
	    if ((comBorders(i).get_L() == process) || (comBorders(i).get_R() == process)) {
			n_borders_loc += 1;
	    }
	}
	
	// Allocate local skeleton
	long nodes_per_border = comBorderNodes.get_n_nodes();
	Skeleton local_skel(n_borders_loc, nodes_per_border, LOCAL);

	// Loop again to generate entries in local_skel
	long border_ix = 0;
	long global_node_ix = 0;
	long local_node_ix = 0;

	for (long i = 0; i < n_borders; ++i) {
	    if ((comBorders(i).get_L() == process) || (comBorders(i).get_R() == process)) {
	    	copy_border_entries(i, local_skel.get_border(border_ix));

		// Copy border nodes and map from global to local node index
		for (long array_ix = 0; array_ix < nodes_per_border; ++array_ix) {
		    global_node_ix = comBorderNodes.get_entry(i, array_ix);
		    bool found = false;
		    for (long j = 0; j < length_l2g; ++j) {
			if (local2global[j] == global_node_ix) {
			    local_node_ix = j;
			    found = true;
			    break;
			}
			if ((j == length_l2g - 1) && (!found)) {
			    std::cerr << "Skeleton.CreateLocal(...): MAPPING FAILED\n" << std::endl;
			    abort();
			}
		    }
		    //write local node ix to local skel
		    local_skel.comBorderNodes.set_entry(border_ix, array_ix, local_node_ix);
		}
		border_ix += 1;
	    }
	}

	// Write local skeleton to current object 
	// Global skeleton locally obsolete
	*this = std::move(local_skel);
    }
} // Namespace Skeleton

