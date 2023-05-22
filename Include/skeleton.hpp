#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "hpc.hpp"

namespace Mesh{
    class Couple {
    private:
    	long index, c1, c2, L, R, color;
    public:       
        void set_entries(long i, long start_node, long end_node, 
                         long left_proc, long right_proc, long couple_color) {
            index = i;
            c1 = start_node;
            c2 = end_node;
            L = left_proc;
            R = right_proc;
            color = couple_color;
        }
        
        void Print();
        
    };
    
    class ICouple {
    private:
    	// all couples are stored as one long list of nodes
        // stuff to refactor as list of lists
    	long n_nodes;	// nodes per couple
    	long n_couples;
        List<long> nodes;
    public:
        ICouple(long n_couples) :
                nodes(n_couples), n_nodes(1), n_couples(n_couples) {}
        
        ICouple(long n_couples, long n_nodes) :
                nodes(n_couples*n_nodes), n_nodes(n_nodes), n_couples(n_couples) {}
                
        long get_n_nodes() {return n_nodes;}
        long get_n_couples() {return n_couples;}
        
        void init_entries(long index) {
            for (long i = 0; i < n_nodes; ++i) {
            	nodes(index*n_nodes + i) = i+1;
            }
        }
                
        void Print();	
    };
    
    class Skeleton {
    private:
    	List<Couple> couples;
        ICouple icouples;
        
        long n_couples, n_icouples;
    public:
        Skeleton(long m, long n) :
                 couples(2*n*m-n-m), icouples(2*n*m-n-m),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
        Skeleton(long m, long n, long refine_factor) :
                 couples(2*n*m-n-m), 
                 icouples(2*n*m-n-m, pow(2, refine_factor) - 1),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
                 
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        long get_n_couples() {return n_couples;}
        long get_n_icouples() {return n_icouples;}
                 
        void Create(Mesh &mesh);
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
