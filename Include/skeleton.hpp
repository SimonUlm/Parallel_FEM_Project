#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "hpc.hpp"

namespace Skeleton{
	
	enum global_or_local { GLOBAL, LOCAL };

    class ComBorder {
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
        
        void copy_entries(ComBorder& border) {
        	border.set_entries(index, c1, c2, L, R, color);
        }
        
        long get_L() {return L;}
        long get_R() {return R;}
        
        void Print();
        
    };
    
    class ComBorderNodes {
    private:
    	// all couples are stored as one long list of nodes
        // stuff to refactor as list of lists
    	long n_nodes;	// nodes per couple
    	long n_borders;
        Util::List<long> nodes;
    public:
        ComBorderNodes(long n_borders) :
                nodes(n_borders), n_nodes(1), n_borders(n_borders) {}
        
        ComBorderNodes(long n_borders, long n_nodes) :
                nodes(n_borders*n_nodes), n_nodes(n_nodes), n_borders(n_borders) {}
                
        long get_n_nodes() {return n_nodes;}
        long get_n_borders() {return n_borders;}
        
        void init_entries(long index) {
            for (long i = 0; i < n_nodes; ++i) {
            	nodes(index*n_nodes + i) = i+1;
            }
        }
        
        void set_entry(long ix_border, long ix_node, long node) {
            nodes(ix_border*n_nodes + ix_node) = node;    
        }
        
        long get_entry(long ix_border, long ix_node) {
            return nodes(ix_border*n_nodes + ix_node);    
        }
                
        void Print();	
    };
    
    
    /* Skeleton Class */
    /*!  \class Skeleton skeleton.hpp "Include/skeleton.hpp"
     *   \brief The skeleton is used for the communciation between the distributed
     *	  processes
     *	 
     *	 The Skeleton class contains a list of the \a ComBorders which contain the
     *	 information of the corresponding processes and the edge nodes. Additionally
     *   Additionally it contains an instance of ComBorderNodes which handles the 
     *   intervening nodes on each comunication border.
     *
     */
    class Skeleton {
    private:
        long n_borders;
        Util::List<ComBorder> comBorders;
        ComBorderNodes comBorderNodes;
        
    public:
    	/*!  
    	 *   Skeleton Constructor for unrefined meshes
    	 */
        Skeleton(long m, long n) :
                 comBorders(2*n*m-n-m), comBorderNodes(2*n*m-n-m),
                 n_borders(2*n*m-n-m) {}
        /*!  
    	 *   Skeleton Constructor for refined meshes
    	 */
        Skeleton(long m, long n, long refine_factor) :
                 comBorders(2*n*m-n-m), 
                 comBorderNodes(2*n*m-n-m, pow(2, refine_factor) - 1),
                 n_borders(2*n*m-n-m) {}
        /*!  
    	 *   Skeleton Constructor for local Skeltons
    	 */
        Skeleton(long n_borders, long nodes_per_border, enum global_or_local usecase) :
				 comBorders(n_borders), comBorderNodes(n_borders * nodes_per_border), 
		         n_borders(n_borders) {assert(usecase == LOCAL);}
                 
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        Skeleton & operator=(Skeleton &&other) = default;
        Skeleton & operator=(const Skeleton &) = delete;
        
        /*!  
    	 *   Return the number of borders of this Skeleton
    	 */
        long get_n_borders() {return n_borders;}
        
        /*!  
    	 *   Function to copy of border to another Skelton (used for creating
    	 *   local Skeletons)
    	 */
        void copy_border_entries(long ix_border, ComBorder& border) {
        	comBorders(ix_border).copy_entries(border);
        }
        
        /*!  
    	 *   Return reference of border (used for copying from global skeleton)
    	 */
        ComBorder& get_border(long ix_border) {
        	return comBorders(ix_border);
        }
        
        /*!  
    	 *   Creating global Skeleton from mesh
    	 */         
        void Create(Mesh::Mesh &mesh);
        
        /*!  
    	 *   Transforms global Skeleton to local Skeleton
    	 */
        void CreateLocal(long process, long* local2global, long length_l2g);
        /*!  
    	 *   Prints Data of Skeleton
    	 */
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
