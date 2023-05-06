#include "hpc.h"


mesh *createMesh(index m, index n)
{
    index ncoord = ((m+1) * (n+1)) * 2;
    index nelem = 2 * m * n;
    index nedges = 3 * m * n + m + n;
    index nbdry = 2 * (m + n);
    index nfixed = m + n + 1;

    mesh* newMesh = mesh_alloc_with_edges(ncoord, nelem, nbdry, nedges, nfixed);

    // Filling in the coordinates
    int k = 0;
    for(index i = 0; i <= m; ++i){
        for(index j = 0; j <= n; ++j){
            double x = i/m; // check if double 
            double y = j/n;
            newMesh->coord[k] = x;
            newMesh->coord[k+1] = y;
            k += 2;
        }

    }

    





}


