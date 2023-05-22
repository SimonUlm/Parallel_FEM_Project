#include <cstdio>

#include "hpc.hpp"

namespace Mesh {
    void Node::Print() {
    	printf("    (%lg,  %lg)\n", x, y);
    }
    
    void Element::Print() {
	printf(" %zu %zu %zu", n1, n2, n3);
        printf(" %zu %zu %zu", m1, m2, m3);
        printf(" %zu", t);
        printf("\n");    
    }
    
    void BoundaryEdge::Print() {
        printf(" %zu %zu", n1, n2);  
        printf(" %zu %zu", m, t);
        printf("\n");
    }
}
