#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "hpc.hpp"

namespace Skeleton{
	
	enum global_or_local {GLOBAL, LOCAL };

    /* ComBorder */
    
    class ComBorder {
    /* 
     * The ComBorder class contains the two crosspoints for this border and the corresponding
     * processes. Additionally a color is saved so border of each process is different
     *
     */
    private:
    	long index; // index of border
        long c1; 	// starting crosspoint
        long c2; 	// end crosspoint
        long L; 	// Left or lower process
        long R; 	// Right or upper process
        long color; // Color of to time communication
    public:
    	// Set entries of ComBorder       
        void set_entries(long i, long start_node, long end_node, 
                         long left_proc, long right_proc, long couple_color) {
            index = i;
            c1 = start_node;
            c2 = end_node;
            L = left_proc;
            R = right_proc;
            color = couple_color;
        }
        
        // Copy entries from another ComBorder
        void copy_entries(ComBorder& border) {
        	border.set_entries(index, c1, c2, L, R, color);
        }
        
        // Getter Methods
        const long get_L() const {return L;}
        const long get_R() const {return R;}
        const long get_c1() const {return c1;}
        const long get_c2() const {return c2;}
        
        // Setter for start/end Node
        void set_c1(long new_c1) {c1 = new_c1;}
        void set_c2(long new_c2) {c2 = new_c2;}
        
        // Print Data of ComBorder
        void Print();
        
    };
    
    /* ComBorderNodes */
    
    class ComBorderNodes {
    /*
     * The ComBorderNodes class consists one Vector of nodes. This list is not seperated into
     * the single borders but is handled as one consectuive Vector. The seperation into each
     * border is handled by the class itself. So from the outside the Borders look seperated
     *
     */
    private:
    	long n_nodes;				// Number of nodes per ComBorder
    	Util::Vector<long> nodes; 	// Consecutive list of all nodes on the ComBorder's
    	
    public:
    	
        /*
    	 * Constructor for unrefined Mesh
         *
         * n_borders: Number of ComBorder's in this Skeleton
         *
    	 */
        ComBorderNodes(long n_borders) :
        	nodes(0), n_nodes(0) {}
        
        /*
    	 * Constructor for refined Mesh
         *
         * n_borders: Number of ComBorder's in this Skeleton
         * n_nodes: Number of nodes per ComBorder (calculated from refinement factor)
         *
    	 */
        ComBorderNodes(long n_borders, long n_nodes) :
        	nodes(n_borders*n_nodes), n_nodes(n_nodes) {}
        
        // Getter Methods       
        long get_n_nodes() { return n_nodes; }
        Util::Vector<long> &get_nodes() { return nodes; }
        
        /*
         * Get node from given border
         *
         * ix_border: index of border
         * ix_node: index of node in border
         *
         */
		long get_entry(long ix_border, long ix_node) {
            return nodes(ix_border*n_nodes + ix_node);    
        }
        
        /*
         * Set node to given border
         *
         * ix_border: index of border
         * ix_node: index of node in border
         *
         */
        void set_entry(long ix_border, long ix_node, long node) {
            nodes(ix_border*n_nodes + ix_node) = node;    
        }
        
        // Print data ComBorderNodes seperated for each border
        void Print();	
    };

    /* Skeleton */
    
    class Skeleton {
    /*
     * The Skeleton is used for communication between the distributed processes.
     * The class contains a vector of ComBorders which contain the information
     * of the corresponding processes and the edge nodes. Additionally it contains
     * an instance of ComBorderNodes which stores all nodes on the communication
     * borders.
     *
     */
    private:
        long n_borders; 						// Number of borders in this Skeleton
        Util::Vector<ComBorder> com_borders;	// Vector of ComBorders
        ComBorderNodes com_border_nodes; 		// List of nodes on the communication borders
        
    public:
    	/*  
    	 * Constructor for unrefined meshes
    	 *
         * m: Number of processes in vertical direction
         * n: Number of processes in horizontal direction
         *
    	 */
        Skeleton(long m, long n) :
                com_borders(2 * n * m - n - m), com_border_nodes(2 * n * m - n - m),
                n_borders(2*n*m-n-m) {}
                
        /*  
    	 * Constructor for refined meshes
         *
         * m: Number of processes in vertical direction
         * n: Number of processes in horizontal direction
         * refine_factor: Number of performed refinements
         *
    	 */
        Skeleton(long m, long n, long refine_factor) :
                com_borders(2 * n * m - n - m),
                com_border_nodes(2 * n * m - n - m, pow(2, refine_factor) - 1),
                n_borders(2*n*m-n-m) {}
                
        /*  
    	 * Constructor for local Skeltons
         *
         * n_borders: Number of ComBorder's
         * nodes_per_border: Number of nodes on ComBorder edge
         *
    	 */
        Skeleton(long n_borders, long nodes_per_border, enum global_or_local usecase) :
                com_borders(n_borders), com_border_nodes(n_borders, nodes_per_border),
                n_borders(n_borders) {assert(usecase == LOCAL);}
               
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        Skeleton & operator=(Skeleton &&other) = default;
        Skeleton & operator=(const Skeleton &) = delete;
        
    	// Creating global Skeleton from Mesh      
        void Create(Mesh::Mesh &mesh);
        
    	//Transforms global Skeleton to local Skeleton
        void CreateLocal(int rank, Mesh::LocalMesh &local_mesh);
        
        #ifdef _MPI
        // Scatter Skeleton between Processes by MPI
        void Scatter(int rank, Mesh::LocalMesh &local_mesh);
        #endif
		
		// Print data of Skeleton
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
