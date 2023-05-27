#include <cstdio>

#include "hpc.hpp"

namespace Mesh {

// Print mesh information
    void Mesh::Print() {

        // Start printing
        printf("\n=========== Print Mesh Data ===========\n");

        // Print mesh sizes
        printf("Mesh size:\n");
        printf("Number of Coordinates: %zu\nNumber of Elements: %zu\nNumber of Boundary-Elements: %zu\n",
               nodes.count, elements.count, boundary.count);

        // Print coordinates
        printf("\nCoordinates (x,y):\n");
        for (auto node : nodes)
            node.Print();

        // Print elements
        printf("\nElements:\n");
        printf("Vertices (n1,n2,n3), Mid. Points (m1,m2,m3), Affiliation\n");
        for (auto element : elements)
            element.Print();

        // Print boundary
        printf("\nBoundary Elements:\n");
        printf("Endpoints (n1, n2), Edge Number (ed1), Type\n");
        for (auto boundary_edge : boundary)
            boundary_edge.Print();

        // If there are fixed_nodes nodes (usually dirichlet nodes) Print them
        if (fixed_nodes.count) {
            printf("\nFixed Nodes:\n");
            for (auto fixed_node : fixed_nodes)
                printf(" %zu\n", fixed_node);
        }

        // Print overall storage requirements
        printf("\nMemory\n");
        printf("Coordinates : %12zu Byte\n", nodes.count * 2 * sizeof(double));
        printf("Elements    : %12zu Byte\n", elements.count * 7 * sizeof(long));
        printf("Boundary    : %12zu Byte\n", boundary.count * 4 * sizeof(long));
        printf("Edge2no     : %12zu Byte\n", edges.count * 2 * sizeof(long));
        long total = nodes.count * 2 * sizeof(double)
                + (7 * elements.count + 4 * boundary.count + edges.count * 2) * sizeof(long);
        printf("Total       : %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
}
