#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "hpc.hpp"

namespace Mesh{
	enum global_or_local { GLOBAL, LOCAL };

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
        
        void copy_entries(Couple& couple) {
        	couple.set_entries(index, c1, c2, L, R, color);
        }
        
        long get_L() {return L;}
        long get_R() {return R;}
        
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
        
        void set_entry(long ix_icouple, long ix_node, long node) {
            nodes(ix_icouple*n_nodes + ix_node) = node;    
        }
        
        long get_entry(long ix_icouple, long ix_node) {
            return nodes(ix_icouple*n_nodes + ix_node);    
        }
                
        void Print();	
    };
    
    class Skeleton {
    private:
        long n_couples, n_icouples;
        List<Couple> couples;
        ICouple icouples;
        
    public:
        Skeleton(long m, long n) :
                 couples(2*n*m-n-m), icouples(2*n*m-n-m),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
        Skeleton(long m, long n, long refine_factor) :
                 couples(2*n*m-n-m), 
                 icouples(2*n*m-n-m, pow(2, refine_factor) - 1),
                 n_couples(2*n*m-n-m), n_icouples(2*n*m-n-m) {}
        Skeleton(long n_couples, long nodes_per_couple, enum global_or_local usecase) :
				 couples(n_couples), icouples(n_couples * nodes_per_couple), 
		         n_couples(n_couples), n_icouples(n_couples) {assert(usecase == LOCAL);}
                 
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        Skeleton & operator=(Skeleton &&other) = default;
        Skeleton & operator=(const Skeleton &) = delete;
        
        long get_n_couples() {return n_couples;}
        long get_n_icouples() {return n_icouples;}
        
        void copy_couple_entries(long ix_couple, Couple& couple) {
        	couples(ix_couple).copy_entries(couple);
        }
        
        Couple& get_couple(long ix_couple) {
        	return couples(ix_couple);
        }
                 
        void Create(Mesh &mesh);
        void CreateLocal(long process, long* local2global, long length_l2g);
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
