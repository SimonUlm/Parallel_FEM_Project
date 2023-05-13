#include "hpc.h"

// Allocate mesh data structure
mesh *mesh_alloc(index ncoord, index nelem, index nbdry) {
    mesh *M = (mesh*) malloc(sizeof(mesh));  /* allocate the mesh struct */
    if (!M) return (NULL);                    /* out of memory */
    M->ncoord = ncoord;
    M->nelem = nelem;
    M->nbdry = nbdry;
    M->nedges = 0;
    M->nfixed = 0;
    M->coord = (double *) malloc(ncoord * 2 * sizeof(double));
    M->elem = (index *) malloc(nelem * 7 * sizeof(index));
    M->bdry = (index *) malloc(nbdry * 4 * sizeof(index));
    M->edge2no = NULL;
    M->fixed = NULL;
    return ((!M->coord || !M->elem || !M->bdry) ? mesh_free(M) : M);
}

// Allocate mesh data structure with edges and fixed_nodes
mesh *mesh_alloc_with_edges(index ncoord, index nelem, index nbdry, index nedges, index nfixed) {
    mesh *M = (mesh *) malloc(sizeof(mesh));  /* allocate the mesh struct */
    if (!M) return (NULL);                    /* out of memory */
    M->ncoord = ncoord;
    M->nelem = nelem;
    M->nbdry = nbdry;
    M->nedges = nedges;
    M->nfixed = nfixed;
    M->coord = (double *) malloc(ncoord * 2 * sizeof(double));
    M->elem = (index *) malloc(nelem * 7 * sizeof(index));
    M->bdry = (index *) malloc(nbdry * 4 * sizeof(index));
    M->edge2no = (index *) malloc(nedges * 2 * sizeof(index));
    M->fixed = (index *) malloc(nfixed * sizeof(index));
    return ((!M->coord || !M->elem || !M->bdry) ? mesh_free(M) : M);
}


// Free a mesh data structure 
mesh *mesh_free(mesh *M) {
    if (!M) return (NULL);     /* do nothing if M already NULL */
    free(M->coord);
    free(M->elem);
    free(M->bdry);
    if (M->edge2no) free(M->edge2no);
    if (M->fixed) free(M->fixed);
    free(M);
    return (NULL);   /* free the mesh struct and return NULL */
}
