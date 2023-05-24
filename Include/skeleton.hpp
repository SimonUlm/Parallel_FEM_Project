#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP
#include <math.h>
#include <cstdio>

#include "hpc.hpp"

namespace Skeleton{
	
	enum global_or_local { GLOBAL, LOCAL };

    /* ComBorder */
    /*!  \class ComBorder skeleton.hpp "Include/skeleton.hpp"
     *   \brief ComBorder manages the crosspoints and processes for each border
     *
     *	 The ComBorder class contains the two crosspoints for this border and the corresponding processes
     *   Additionally a color is saved so each border of each process is different
     */
    class ComBorder {
    private:
    	long index; /*!< index of border */
        long c1; /*!< starting crosspoint */
        long c2; /*!< end crosspoint */
        long L; /*!< Left or lower process */
        long R; /*!< Right or upper process */
        long color; /*!< Color of to time communication */
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
        
        /*!
    	 *   \brief Prints Data of ComBorder
    	 */
        void copy_entries(ComBorder& border) {
        	border.set_entries(index, c1, c2, L, R, color);
        }
        
        long get_L() {return L;}
        long get_R() {return R;}
        
        /*!
    	 *   Prints Data of ComBorder
    	 */
        void Print();
        
    };
    
    /* ComBorderNodes */
    /*!  \class ComBorderNodes skeleton.hpp "Include/skeleton.hpp"
     *   \brief ComBorderNodes manages all corresponding nodes on all communciation borders of the skeleton
     *
     *	 The ComBorderNodes class consists a List of nodes. This list is not seperated into single ComBorder's but is
     *   handled as a consectuive list. The seperation into each border is handled by the setter and getter of this class
     *
     */
    class ComBorderNodes {
    private:
    	// all couples are stored as one long list of nodes
        // stuff to refactor as list of lists
    	long n_nodes;	/*!< Number of nodes per ComBorder */
    	long n_borders; /*!< Number of ComBorder's */
        
    public:
    	Util::List<long> nodes; /*!< Consectuive list of all nodes on the ComBorder's */
        /*!
    	 *   Constructor for unrefined Mesh (Check if still working with 0??)
         *
         *   \param n_borders Number of ComBorder's in this Skeleton
    	 */
        ComBorderNodes(long n_borders) :
                nodes(0), n_nodes(0), n_borders(n_borders) {}
        
        /*!
    	 *   Constructor for refined Mesh
         *
         *   \param n_borders Number of ComBorder's in this Skeleton
         *   \param n_nodes Number of nodes per ComBorder (calculated from refinement factor)
    	 */
        ComBorderNodes(long n_borders, long n_nodes) :
                nodes(n_borders*n_nodes), n_nodes(n_nodes), n_borders(n_borders) {}
                
        long get_n_nodes() {return n_nodes;}
        long get_n_borders() {return n_borders;}

        /*!
    	 *   Initialize nodes with ascending numbers per border
    	 */
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

        /*!
    	 *   \brief Prints Data of ComBorderNodes
         *
         *   The nodes for each border are printed into a new line
    	 */
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
        long n_borders; /*!< Number of borders in this Skeleton */
        Util::List<ComBorder> comBorders; /*!< List of ComBorder */
        ComBorderNodes comBorderNodes; /*!< ComBorderNodes consist a list of nodes corresponding to each border */
        
    public:
    	/*!  
    	 *   Skeleton Constructor for unrefined meshes
         *   \param m Number of processes in vertical direction
         *   \param n Number of processes in horizontal direction
    	 */
        Skeleton(long m, long n) :
                 comBorders(2*n*m-n-m), comBorderNodes(2*n*m-n-m),
                 n_borders(2*n*m-n-m) {}
        /*!  
    	 *   Skeleton Constructor for refined meshes
         *
         *   \param m Number of processes in vertical direction
         *   \param n Number of processes in horizontal direction
         *   \param refine_factor Number of performed refinements
    	 */
        Skeleton(long m, long n, long refine_factor) :
                 comBorders(2*n*m-n-m), 
                 comBorderNodes(2*n*m-n-m, pow(2, refine_factor) - 1),
                 n_borders(2*n*m-n-m) {}
        /*!  
    	 *   Skeleton Constructor for local Skeltons
         *
         *   \param n_borders Number of ComBorder's
         *   \param nodes_per_border Number of nodes on ComBorder edge
    	 */
        Skeleton(long n_borders, long nodes_per_border, enum global_or_local usecase) :
				 comBorders(n_borders), comBorderNodes(n_borders, nodes_per_border),
		         n_borders(n_borders) {assert(usecase == LOCAL);}
                 
        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        Skeleton & operator=(Skeleton &&other) = default;
        Skeleton & operator=(const Skeleton &) = delete;
        
        /*!  
    	 *   \brief Return the number of borders of this Skeleton
    	 */
        long get_n_borders() {return n_borders;}
        
        /*!  
    	 *   \brief Function to copy of ComBorder entries
         *
         *   This function is used in CreateLocal to copy values of the ComBorder's to the tempoary local
         *   Skeleton
    	 */
        void copy_border_entries(long ix_border, ComBorder& border) {
        	comBorders(ix_border).copy_entries(border);
        }
        
        /*!  
    	 *   \brief Get a reference of a ComBorder
         *
         *   This function is used to sent a ComBorder from the tempoary local Skeleton to the global one
         *   to copy entries from the global Skeleton to the local Skeleton
         *
         *   \param ix_border Index of ComBorder
         *   \return Refererence of ComBorder with given index
    	 */
        ComBorder& get_border(long ix_border) {
        	return comBorders(ix_border);
        }
        
        void set_border_node(long ix_border, long ix_node, long node) {
        	comBorderNodes.set_entry(ix_border, ix_node, node);
        }
        /*!  
    	 *   \brief Creating global Skeleton from Mesh
         *
         *   Creates the global Skeleton from given Mesh. With the number of processes of the Mesh
         *   the Crosspoints and Processes of the ComBorder's are calculated. With the refinement
         *   factor the ComBorderNodes are calcualted and saved.
         *
         *   \param mesh Mesh for which the Skeleton is created
         *
    	 */         
        void Create(Mesh::Mesh &mesh);
        
        /*!  
    	 *   \brief Transforms global Skeleton to local Skeleton
         *
         *   Creates a tempory Skeleton copies the ComBorder's which correspond with the local process
         *   to the tempory Skeleton. At the end the old global Skeleton is overwritten with the local one
         *
         *   \param process Number of Local Process
         *   \param local2global
         *   \param langth_l2g
    	 */
        void CreateLocal(int rank, Mesh::LocalMesh &local_mesh);
        
        #ifdef _MPI
        void Scatter(int rank, Mesh::LocalMesh &local_mesh);
        #endif

        /*!
    	 *   Prints Data of Skeleton
    	 */
        void Print();
    };
}

#endif //HPC2_SKELETON_HPP
