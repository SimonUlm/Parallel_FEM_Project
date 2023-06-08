#include <cstdio>

#include "hpc.hpp"

namespace Mesh {
    // Print x and y coordinates of the Node
    void Node::Print() {
        printf("    (%lg,  %lg)\n", x, y);
    }

    // Print all Nodes, all Edges and the Affiliation belonging to the Element
    void Element::Print() {
        printf(" %zu %zu %zu", n1, n2, n3);
        printf(" %zu %zu %zu", m1, m2, m3);
        printf(" %zu", t);
        printf("\n");
    }

    // Print Nodes, Edge number and Type of the Boundary Edge
    void BoundaryEdge::Print() {
        printf(" %zu %zu", n1, n2);
        printf(" %zu %zu", m, t);
        printf("\n");
    }
}
