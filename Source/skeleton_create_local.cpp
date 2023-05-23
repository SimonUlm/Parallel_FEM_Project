#include "hpc.hpp"

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
	Skeleton local_skel(n_borders_loc, comBorderNodes.get_n_nodes(), LOCAL);

	// Loop again to generate entries in local_skel
	long border_index = 0;
	for (long i = 0; i < n_borders; ++i) {
	    if ((comBorders(i).get_L() == process) || (comBorders(i).get_R() == process)) {
	    	copy_border_entries(i, local_skel.get_border(border_index));

			border_index += 1;
	    }
	}

	// Write local skeleton to current object 
	// Global skeleton locally obsolete
	*this = std::move(local_skel);
    }
} // Namespace Skeleton

