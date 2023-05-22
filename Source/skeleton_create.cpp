#include "hpc.hpp"

namespace Mesh {

    void Skeleton::Create(Mesh &mesh) {  
        long couple_index = 0;
        long m = mesh.get_m();
        long n = mesh.get_n();

        long nof_h_edges = (m + 1) * n; // number of horizontal edges
        long nof_v_edges = m * (n + 1); // number of vertical edges
        long nodes_per_row = n + 1;
        
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j < n; ++j) {
                if (j < n-1) {
            	    // initialize right vertical couple
            	    long c1 = nodes_per_row * i + j + 1;
            	    long c2 = nodes_per_row * (i + 1) + (j + 1);
            	    long l_proc = i * n + j;
            	    long r_proc = i * n + j + 1;
            	    long color = j%2;
            	    couples(couple_index).set_entries(
            	    	couple_index, c1, c2, l_proc, r_proc, color
            	    );
            	    icouples.init_entries(couple_index);     	    
            	    couple_index++; 	    
            	}
            	if (i < m-1) {
            	    // initialize top horizontal couple
            	    long c1 = nodes_per_row * (i + 1) + (j + 1);
            	    long c2 = nodes_per_row * (i + 1) + (j + 1) - 1;
            	    long l_proc = i * n + j;
            	    long r_proc = (i + 1) * n + j;
            	    long color = i%2 + 2;
            	    couples(couple_index).set_entries(
            	    	couple_index, c1, c2, l_proc, r_proc, color
            	    );
            	    icouples.init_entries(couple_index);   	    
            	    couple_index++; 	    
            	}
            
            }
        }
    }
}
