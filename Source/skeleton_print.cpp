#include <cstdio>

#include "skeleton.hpp"

namespace Mesh {

// Print mesh information
    void Skeleton::Print() {
        long brief = 0;
        // Start printing
        printf("\n=========== Print Skeleton Data ===========\n");
        
        long n_nodes = icouples.n_nodes;
        // Print mesh sizes
        printf("Skeleton size:\n");
        printf("Number of Couples: %zu\nNumber of ICouples: %zu\nNumber of Nodes per ICouple: %zu\n",
               n_couples, n_icouples, n_nodes);
               
        // Print Couples
        printf("\nCouples:\n");
        printf("Index, Start_Node, End_Node, L-Proc, R-Proc, Color\n");
        for (long i = 0; i < n_couples; ++i) {
            couples(i).Print();
        }
        
        // Print ICouples
        printf("\nICouples:\n");
        printf("Connecting Nodes\n");
        icouples.Print();

        // Print overall storage requirements
        printf("\nMemory\n");
        printf("Couples     : %12zu Byte\n", n_couples * 6 * sizeof(long));
        printf("ICouple     : %12zu Byte\n", n_icouples * n_nodes * sizeof(long));
        long total = n_couples * 6 * sizeof(long) + n_icouples * n_nodes * sizeof(long);
        printf("Total       : %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
}
