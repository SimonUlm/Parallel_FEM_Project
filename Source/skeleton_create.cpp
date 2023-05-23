#include <math.h>
#include "hpc.hpp"

namespace Mesh {
    long mn_r(long m, long r);

    void Skeleton::Create(Mesh &mesh) {  
        long couple_index = 0;
        long m = mesh.get_m();
        long m_n = m + 1;
        long n = mesh.get_n();
        long n_n = n + 1;

        long nof_h_edges = (m + 1) * n; // number of horizontal edges
        long nof_v_edges = m * (n + 1); // number of vertical edges
        long nodes_per_row = n + 1;
        long refine = mesh.get_refine_factor();
        
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j < n; ++j) {
            	if (j != 0) {
            	    // initialize left vertical couple
            	    long c1 = i * nodes_per_row + j;		// lower left Node
            	    long c2 = nodes_per_row * (i + 1) + j;	// upper left Node
            	    long l_proc = i * n + (j-1);		// left process
            	    long r_proc = i * n + j;			// this process
            	    long color = (j-1)%2;			// color
            	    
            	    couples(couple_index).set_entries(
            	    	couple_index, c1, c2, l_proc, r_proc, color
            	    );
            	    
            	    long index = 0;
            	    for (long r = refine; r > 0; --r) {
            	        long node = mn_r(m_n, r-1) * mn_r(n_n, r-1) + (mn_r(n_n, r-1) -1) * m_n;
            	        node += pow(2, r-1) * j; 		// accounting column
            	        node += pow(2, r-1) * n_n * i;		// accounting row
            	        for (long k = 0; k < pow(2, r-1); ++k) {
			    icouples.set_entry(couple_index, index, node);
            	            node += 1;
            	            index += 1;
            	        }    
            	    }   
            	    couple_index++;	
            	}
            	
            	if (i != 0) {
            	    // initialize low horizontal couple
            	    long c1 = nodes_per_row * i + j + 1;	// lower right node
            	    long c2 = nodes_per_row * i + j;		// lower left node
            	    long l_proc = (i - 1) * n + j;		// lower process
            	    long r_proc = i * n + j;			// this process as upper process
            	    long color = (i-1)%2 + 2;			// color
            	    
            	    couples(couple_index).set_entries(
            	    	couple_index, c1, c2, l_proc, r_proc, color
            	    );
            	    
            	    long index = 0;
            	    for (long r = refine; r > 0; --r) {
            	        long node = mn_r(m_n, r-1) * mn_r(n_n, r-1);
            	        node += pow(2, r-1) * j; 		// accounting column
            	        node += pow(2, r-1) * n * i;		// accounting row
            	        for (long k = 0; k < pow(2, r-1); ++k) {
			    icouples.set_entry(couple_index, index, node);
            	            node += 1;
            	            index += 1;
            	        }    
            	    }   
            	    couple_index++;	
            	}
            }
        }
    }
    
    long mn_r(long m, long r) {
    	// m = number of nodes in one dimension
    	return (m - 1) * pow(2, r) + 1;
    }
}
