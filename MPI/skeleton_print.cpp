#include <cstdio>

#include "hpc.hpp"

namespace Skeleton {

    /*
     * Prints Data of Skeleton
     */
    void Skeleton::Print() {
        long brief = 0;
        // Start printing
        printf("\n=========== Print Skeleton Data ===========\n");

        long n_nodes = com_border_nodes_.n_nodes();
        // Print mesh sizes
        printf("Skeleton size:\n");
        printf("Number of ComBorders: %zu\nNumber of Nodes per ComBorder: %zu\n",
               n_borders_, n_nodes);

        // Print Couples
        printf("\nComBorders:\n");
        printf("Index, Start_Node, End_Node, L-Proc, R-Proc, Color\n");
        for (long i = 0; i < n_borders_; ++i) {
            com_borders_(i).Print();
        }

        // Print ICouples
        printf("\nComBorderNodes:\n");
        printf("Connecting Nodes\n");
        com_border_nodes_.Print();

        // Print overall storage requirements
        printf("\nMemory\n");
        printf("ComBorders  	  : %12zu Byte\n", n_borders_ * 6 * sizeof(long));
        printf("ComBorderNodes    : %12zu Byte\n", n_borders_ * n_nodes * sizeof(long));
        long total = n_borders_ * 6 * sizeof(long) + n_borders_ * n_nodes * sizeof(long);
        printf("Total       : %12.6g MByte\n", (double) total / 1024. / 1024.);
    }

    /*
     * Prints data ComBorder
     */
    void ComBorder::Print() {
        printf(" %zu %zu %zu", index_, c1_, c2_);
        printf(" %zu %zu %zu", L_, R_, color_);
        printf("\n");
    }

    /*
     * Prints ComBorder Nodes seperated for each ComBorder
     */
    void ComBorderNodes::Print() {
        for (long i = 0; i < nodes_.count() / n_nodes_; ++i) {
            for (long j = 0; j < n_nodes_; ++j) {
                printf(" %zu", nodes_(i * n_nodes_ + j));
            }
            printf("\n");
        }
    }
}
