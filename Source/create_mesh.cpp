#include "hpc.h"


mesh *createMesh(index m, index n)
{
    index ncoord = ((m+1) * (n+1)) * 2;
    index nelem = 2 * m * n;
    index nedges = 3 * m * n + m + n;
    index nbdry = 2 * (m + n);
    index nfixed = m + n + 1;

    mesh* newMesh = mesh_alloc_with_edges(ncoord, nelem, nbdry, nedges, nfixed);

    index indexCoord = 0;

    index hEdges = (m + 1) * n; // number of horizontal edges
    index vEdges = m * (n + 1);
    index nodesPerRow = n + 1;

    // Iterating over rectangles (2 elements at a time)
    for (index i = 0; i <= m; ++i) {
        for (index j = 0; i <= n; ++j) {
            // Calculating node coordinates
            double x = (1.0*i)/m;
            double y = (1.0*j)/n;
            newMesh->coord[indexCoord] = x;
            newMesh->coord[indexCoord+1] = y;
            indexCoord += 2;

            // Calculating the four nodes of each rectangle
            index n1 = i * nodesPerRow + j;                 // ll node
            index n2 = i * nodesPerRow + (j + 1);           // lr node
            index n3 = (i + 1) * nodesPerRow + (j + 1);     // ur node
            index n4 = (i + 1) * nodesPerRow + j;           // ul node

            // Calculating the five edges of each rectangle
            index m1 = i * n + j;                           // lh
            index m2 = (i + 1) * n + j;                     // uh
            
            index m3 = hEdges + i * nodesPerRow + j;        // lv
            index m4 = hEdges + i * nodesPerRow + (j + 1);  // rv    

            index m5 = hEdges + vEdges + i * n + j;         // diag

            index t1 = 0;

            index indexElem = (i * n + j) * 2;

            writeElement(newMesh, indexElem, n1, n2, n3, m1, m4, m5, t1);
            writeElement(newMesh, indexElem+1, n1, n3, n4, m5, m2, m3, t1); 

            // Calculating boundaries
            if (i == 0) {
                index indexBdry = j;
                writeBdry(newMesh, indexBdry, n1, n2, m1, t1);
                writeEdgeNo(newMesh, m1, n1, n2);
                writeFixed(newMesh, n1, n2, m1, t1);
            }
            if (i == m) {
                index indexBdry = n + 2 * m + j;
                writeBdry(newMesh, indexBdry, n3, n4, m2, t1);
            }
            if (j == 0) {
                index indexBdry = n + i * 2;
                writeBdry(newMesh, indexBdry, n1, n4, m3, t1);
                writeEdgeNo(newMesh, m3, n1, n4);
                writeFixed(newMesh, n1, n4, m3, t1);
            }
            if (j == n) {
                index indexBdry = n + i * 2 + 1;
                writeBdry(newMesh, indexBdry, n2, n3, m4, t1);
            }

            writeEdgeNo(newMesh, m2, n3, n4);
            writeEdgeNo(newMesh, m4, n2, n3);
            writeEdgeNo(newMesh, m5, n1, n3);
        }
    }
}


void writeElement(mesh* wMesh, index elem, 
                  index n1, index n2, index n3, 
                  index m1, index m2, index m3, 
                  index t1)
{
    wMesh->elem[elem * 7]     = n1;
    wMesh->elem[elem * 7 + 1] = n2;
    wMesh->elem[elem * 7 + 2] = n3;
    wMesh->elem[elem * 7 + 3] = m1;
    wMesh->elem[elem * 7 + 4] = m2;
    wMesh->elem[elem * 7 + 5] = m3;
    wMesh->elem[elem * 7 + 6] = t1;
}

void writeBdry(mesh* wMesh, index bdry, 
               index n1, index n2, 
               index m1, index t1)
{
    wMesh->bdry[bdry * 4]     = n1;
    wMesh->bdry[bdry * 4 + 1] = n2;
    wMesh->bdry[bdry * 4 + 2] = m1;
    wMesh->bdry[bdry * 4 + 3] = t1;
}

void writeEdgeNo(mesh* wMesh, index edge, index n1, index n2)
{
    wMesh->edge2no[edge * 2]     = n1;
    wMesh->edge2no[edge * 2 + 1] = n2;
}

void writeFixed(mesh* wMesh, 
                index n1, index n2, 
                index m1, index t1)
{
    static index fixed = 0;
    wMesh->fixed[fixed * 4] = n1;
    wMesh->fixed[fixed * 4 + 1] = n2;
    wMesh->fixed[fixed * 4 + 2] = m1;
    wMesh->fixed[fixed * 4 + 3] = t1;
    fixed++;

}


