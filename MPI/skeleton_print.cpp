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
        
        long n_nodes = com_border_nodes.get_n_nodes();
        // Print mesh sizes
        printf("Skeleton size:\n");
        printf("Number of ComBorders: %zu\nNumber of Nodes per ComBorder: %zu\n",
               n_borders, n_nodes);
               
        // Print Couples
        printf("\nComBorders:\n");
        printf("Index, Start_Node, End_Node, L-Proc, R-Proc, Color\n");
        for (long i = 0; i < n_borders; ++i) {
            com_borders(i).Print();
        }
        
        // Print ICouples
        printf("\nComBorderNodes:\n");
        printf("Connecting Nodes\n");
        com_border_nodes.Print();

        // Print overall storage requirements
        printf("\nMemory\n");
        printf("ComBorders  	  : %12zu Byte\n", n_borders * 6 * sizeof(long));
        printf("ComBorderNodes    : %12zu Byte\n", n_borders * n_nodes * sizeof(long));
        long total = n_borders * 6 * sizeof(long) + n_borders * n_nodes * sizeof(long);
        printf("Total       : %12.6g MByte\n", (double) total / 1024. / 1024.);
    }
    
    /*
     * Prints Data ComBorder
     */
    void ComBorder::Print() {
    	printf(" %zu %zu %zu", index, c1, c2);
        printf(" %zu %zu %zu", L, R, color);
        printf("\n");
    }
    
    /*
     * Prints ComBorder Nodes seperated for each ComBorder
     */
    void ComBorderNodes::Print() {
    	for (long i = 0; i < nodes.count() / n_nodes; ++i) {
            for(long j = 0; j < n_nodes; ++j) {
            	printf(" %zu", nodes(i*n_nodes+j));
            }
            printf("\n");
        }
    } 
} // Namespace Skeleton
