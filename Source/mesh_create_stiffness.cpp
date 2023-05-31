#include <cassert>
#include <algorithm>
#include <utility>

#include "hpc.hpp"

namespace Mesh {

    Util::SedMatrix Mesh::CreateStiffness() {
        // Missing asserts?
        
        long col = 0;
        long row = 0;

        static int ai[3] = {0,0,1}, aj[3] = {1, 2, 2};
        
        long nElem = elements.count();
        long n = nodes.count();
        
        long ind[3];
        
        long *tmp = new long[n](); // general workspace
        long *ptr_ind = new long[n + 1 + 3 * nElem]();
        
        
        // Count entries per column
        // tmp counts entries per column
        for (long k = 0; k < nElem; ++k) {
        	ind[0] = elements(k).n1;
        	ind[1] = elements(k).n2;
        	ind[2] = elements(k).n3;
        	
        	for (j = 0 ; j < 3 ; j++) {
        		// Get col of entry
        		col = std::min(ind[ai[j]], ind[aj[j]])
            	tmp[col] += 1;
        	}
        	
        }
        
        // Column pointers
        // tmp is set to column offsets in data array of matrix
        long nz = n + 1;
        for (long i = 0; i < n; ++i) {
        	ptr_ind[i] = nz;
        	nz += tmp[i];
        	tmp[i] = ptr_ind[i];        
        }
        ptr_ind[n] = nz;
        
        // Insert indices
        // tmp contains column offsets which are individually increment
        //     to fill in all entries of column on an other position
        for (long k = 0; k < nElem; ++k) {
        	ind[0] = elements(k).n1;
        	ind[1] = elements(k).n2;
        	ind[2] = elements(k).n3;
        	for (j = 0 ; j < 3 ; j++){
        		// Get Col of entry
		        col = std::min(ind[ai[j]], ind[aj[j]]);
		        // Get Row of entry
		        row = std::max(ind[ai[j]], ind[aj[j]]);
		        
		        ptr_ind[tmp[col]]= row;
		        
		        // Write at next position for another row in this column
		        tmp[col] += 1;
        	}
        }
        
        
        // remove duplicate entries
        // tmp checks if row i is yet seen
        long *entries_per_col = new long[n]();
        for (long i = 0 ; i < n ; i++) tmp[i] = n ;  
        
        nz = ptr_ind[0];
        
        for (long j = 0; j < n; ++j) {
        	long q = nz; // column j will start at position q in array
        	
        	// Iterate through elements of this column
        	for (p = ptr_ind[j]; p < ptr_ind[j+1]; ++p) {
        		long i = ptr_ind[p];
        		if !(w[i] >= q) {
        			// none duplicate
        			w[i] = nz;
        			ptr_ind[nz] = i;
        			nz += 1;
        		}    	
        	}
        	ptr_ind[j] = q;
        }
        ptr_ind[n] = nz;
        
        Util::SedMatrix mtrx(n, nz);
        
        for (long j = 0; j < nz + 1; ++j) {
        	mtrx.set_ptr(j, ptr_ind[j]);
        }
        
        return std::move(mtrx);
    }
    
    
}



