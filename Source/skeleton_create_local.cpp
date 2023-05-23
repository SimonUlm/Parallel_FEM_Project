#include "hpc.hpp"

namespace Mesh {
    void Skeleton::CreateLocal(long process, long* local2global, 
			       long length_l2g) {
	// Determine size of local skeleton
	long n_couples_loc = 0;

	for (long i = 0; i < n_couples; ++i) {
	    if ((couples(i).get_L() == process) || (couples(i).get_R() == process)) {
			n_couples_loc += 1;
	    }
	}
	
	// Allocate local skeleton
	Skeleton local_skel(n_couples_loc, icouples.get_n_nodes(), LOCAL);

	// Loop again to generate entries in local_skel
	long couple_index = 0;
	for (long i = 0; i < n_couples; ++i) {
	    if ((couples(i).get_L() == process) || (couples(i).get_R() == process)) {
	    	copy_couple_entries(i, local_skel.get_couple(couple_index));

			couple_index += 1;
	    }
	}

	// Write local skeleton to current object 
	// Global skeleton locally obsolete
	*this = std::move(local_skel);
    }
} // Namespace Mesh

