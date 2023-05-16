#include "mesh.hpp"

namespace Mesh {

// Print mesh information
    void RectangularMesh::print() {
        long brief = 0;

        // Aux variables
        long j, k, ncoord, nelem, nbdry, nfixed, *Elem, *Bdry, *Fixed, total;
        double *Coord;

        // Assign mesh attributes locally
        ncoord = nodes.count;
        nelem = elements.count;
        nbdry = boundary.count;
        nfixed = fixed_nodes.count;
        Coord = &nodes(0).x;
        Elem = &elements(0).n1;
        Bdry = &boundary(0).n1;
        Fixed = &fixed_nodes(0);

        // Start printing
        printf("\n=========== Print Mesh Data ===========\n");

        // Print mesh sizes
        printf("Mesh size:\n");
        printf("Number of Coordinates: %zu\nNumber of Elements: %zu\nNumber of Boundary-Elements: %zu\n",
               ncoord, nelem, nbdry);

        // Print coordinates
        printf("\nCoordinates (x,y):\n");
        for (j = 0; j < ncoord; j++) {
            printf("    (%lg,  %lg)\n", Coord[2 * j], Coord[2 * j + 1]);
            if (brief && j > 10) {
                printf("  ...\n");
                break;
            }
        }

        // Print elements
        printf("\nElements:\n");
        printf("Vertices (n1,n2,n3), Mid. Points (m1,m2,m3), Affiliation\n");
        for (j = 0; j < nelem; j++) {
            for (k = 0; k < 7; k++) {
                printf(" %zu", Elem[7 * j + k]);
            }
            printf("\n");
            if (brief && j > 10) {
                printf("  ...\n");
                break;
            }
        }

        // Print boundary
        printf("\nBoundary Elements:\n");
        printf("Endpoints (n1, n2), Edge Number (ed1), Type\n");
        for (j = 0; j < nbdry; j++) {
            for (k = 0; k < 4; k++) {
                printf(" %zu", Bdry[4 * j + k]);
            }
            printf("\n");
            if (brief && j > 10) {
                printf("  ...\n");
                break;
            }
        }

        // If there are fixed_nodes nodes (usually dirichlet nodes) print them
        if (nfixed) {
            printf("\nFixed Nodes:\n");
            for (j = 0; j < nfixed; j++) {

                printf(" %zu\n", Fixed[j]);

                if (brief && j > 10) {
                    printf("  ...\n");
                    break;
                }
            }
        }

        // Print overall storage requirements
        printf("\nMemory\n");
        printf("Coordinates : %12zu Byte\n", ncoord * 2 * sizeof(double));
        printf("Elements    : %12zu Byte\n", nelem * 7 * sizeof(long));
        printf("Boundary    : %12zu Byte\n", nbdry * 4 * sizeof(long));
        printf("Edge2no     : %12zu Byte\n", edges.count * 2 * sizeof(long));
        total = ncoord * 2 * sizeof(double)
                + (7 * nelem + 4 * nbdry + edges.count * 2) * sizeof(long);
        printf("Total       : %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
}
