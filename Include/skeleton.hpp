#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "mesh.hpp"
#include "mesh_objects.hpp"
#include "mesh_list.hpp"

namespace Mesh{
    class Couple {
    public:
        long index, c1, c2, L, R, color;
               
        void set_entries(long i, long start_node, long end_node, 
                         long left_proc, long right_proc, long couple_color) {
            index = i;
            c1 = start_node;
            c2 = end_node;
            L = left_proc;
            R = right_proc;
            color = couple_color;
        }
        
        void Print() {
            printf(" %zu", index);
            printf(" %zu", c1);
            printf(" %zu", c2);
            printf(" %zu", L);
            printf(" %zu", R);
            printf(" %zu \n", color);
        }
        
    };
    
    class ICouple {
    public:
        // all couples are stored as one long list of nodes
        // stuff to refactor as list of lists
    	long n_nodes;	// nodes per couple
    	long n_couples;
        List<Node> nodes;
        
        ICouple(long n_couples) :
                nodes(n_couples), n_nodes(1), n_couples(n_couples) {}
        
        ICouple(long n_couples, long n_nodes) :
                nodes(n_couples*n_nodes), n_nodes(n_nodes), n_couples(n_couples) {}
                
        void Print() {
            for (long i = 0; i < n_couples; ++i) {
                for(long j = 0; j < n_nodes; ++j) {
                    printf(" %zu", nodes(i*n_nodes+j));
                }
                printf("\n");
            } 
            
        }	
    };
    
    class Skeleton {
    public:
        List<Couple> couples;
        ICouple icouples;
        
        long n_couples, n_icouples;
        
        
        Skeleton(long m, long n) :
                 couples(2*n*m-n-m), icouples(2*n*m-n-m),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
        Skeleton(long m, long n, long refine_factor) :
                 couples(2*n*m-n-m), 
                 icouples(2*n*m-n-m, pow(2, refine_factor) - 1),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
                 
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
                 
        void Create(Mesh &mesh);
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
