#include <cassert>
#include <algorithm>
#include <utility>

#include "hpc.hpp"
#include "mesh_objects.hpp"

namespace Mesh {
    double kappa(Node &n, long typ) {
		return ( 1.0 );
	}
	
    
    void stima_laplace(Node &n1, Node &n2, Node &n3, 
    				   long  typ, double dx[3], double ax[3]) {

    	Node n32 = n3 - n2;
    	Node n13 = n1 - n3;
    	Node n21 = n2 - n1;
    	
    	double fac = (kappa(n1,typ) + kappa(n2,typ) + kappa(n3,typ) ) / 
		      (6.0*(n32.x*n13.y-n32.y*n13.x));
		      
		// Get super diagonal entries
		ax[0] = fac * n13.dot(n32);
		ax[1] = fac * n21.dot(n32);
		ax[2] = fac * n21.dot(n13);

		// Get diagonal entries
		dx[0] = fac * n32.dot(n32);
		dx[1] = fac * n13.dot(n13);
		dx[2] = fac * n21.dot(n21);			    					   
    } 

    Util::SedMatrix Mesh::CreateStiffness() {
        // Missing asserts?
        
        long col = 0;
        long row = 0;
		
		// Used to compare nodes of element to decide in which row of
		// column a entry is needed
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
        	
        	for (long j = 0 ; j < 3 ; j++) {
        		// Get col of entry
        		col = std::min(ind[ai[j]], ind[aj[j]]);
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
        	for (long j = 0 ; j < 3 ; j++){
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
        	for (long p = ptr_ind[j]; p < ptr_ind[j+1]; ++p) {
        		long i = ptr_ind[p];
        		if (!(tmp[i] >= q)) {
        			// none duplicate
        			tmp[i] = nz;
        			ptr_ind[nz] = i;
        			nz += 1;
        		}    	
        	}
        	ptr_ind[j] = q;
        }
        ptr_ind[n] = nz;
        
        
        // Create SedMatrix and initialize non zero pattern
        Util::SedMatrix mtrx(n, nz);
        
        for (long j = 0; j < nz + 1; ++j) {
        	mtrx.set_ptr(j, ptr_ind[j]);
        }
        
        // Element stiffness diagonal and super diagonal
    	double dx[3], ax[3];
    	
        // Get node number per element
        for (long k = 0; k < nElem; ++k) {
        	ind[0] = elements(k).n1;
        	ind[1] = elements(k).n2;
        	ind[2] = elements(k).n3;
        	
        	stima_laplace(nodes(ind[0]), nodes(ind[1]), nodes(ind[2]), elements(k).t, dx, ax);
        	
        	// Add local stiffness diag to global
		    for (long j = 0 ; j < 3 ; ++j){
		        mtrx.add_val(ind[j], ind[j], dx[j]);
		    }
		    
		    
		    for (long j = 0 ; j < 3 ; j++){
        		// Get Col of entry
		        col = std::min(ind[ai[j]], ind[aj[j]]);
		        // Get Row of entry
		        row = std::max(ind[ai[j]], ind[aj[j]]);
		        
		        mtrx.add_val(row, col, ax[j]);
		    }
        }
        return std::move(mtrx);
    }  
}



