#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP

#ifdef _MPI // only needed for parallel algorithms
#include <cstdio>
#include <math.h>
#include <mpi.h>

#include "hpc.hpp"

namespace Skeleton{

    constexpr long kNumberOfColors = 4;     // number of possible colors
	enum global_or_local { GLOBAL, LOCAL }; // Enum for local skeleton creation


    /* ComBorder */

    class ComBorder {
    /*
     * The ComBorder class contains the two crosspoints for this border and the
     * corresponding processes. Additionally, a color is saved so each border of
     * process has a different color
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
        
        // General getter methods
        const long get_L() const {return L;}
        const long get_R() const {return R;}
        const long get_c1() const {return c1;}
        const long get_c2() const {return c2;}
        const long get_color() const { return color; }

        // Setter for start/end Node
        void set_c1(long new_c1) {c1 = new_c1;}
        void set_c2(long new_c2) {c2 = new_c2;}
        void set_index(long new_ix) {index = new_ix;}

        // Print Data of ComBorder
        void Print();
        
    };
    

    /* ComBorderNodes */

    class ComBorderNodes {
    /*
     * The ComBorderNodes class consists one Vector of nodes. This vector is not
     * seperated into the single borders but is handled as one consecutive Vector.
     * The separation into each border is handled by the class itself. So from the
     * outside the Borders look seperated
     *
     */
    private:
    	long n_nodes;				// Number of nodes per ComBorder
    	Util::Vector<long> nodes; 	// Consecutive list of nodes on the ComBorder's

    public:

        /*
    	 * Constructor for unrefined Mesh
         *
         * n_borders: Number of ComBorder's in this Skeleton
         *
    	 */
        ComBorderNodes(long n_borders) :
        	           nodes(0),
        	           n_nodes(0) {}
        
        /*
    	 * Constructor for refined Mesh
         *
         * n_borders: Number of ComBorder's in this Skeleton
         * n_nodes: Number of nodes per ComBorder (calculated from refinement factor)
         *
    	 */
        ComBorderNodes(long n_borders, long n_nodes) :
        	           nodes(n_borders*n_nodes),
        	           n_nodes(n_nodes) {}

        // General getter methods
        const long get_n_nodes() const { return n_nodes; }
        const Util::Vector<long> &get_nodes() const { return nodes; }

        /*
         * Get/Set node for given border
         *
         * ix_border: index of border
         * ix_node: index of node in border
         *
         */
		const long get_entry(long ix_border, long ix_node) const {
            return nodes(ix_border*n_nodes + ix_node);
        }

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
        long n_borders; 					    // Number of borders in this Skeleton
        Util::Vector<ComBorder> com_borders;    // Vector of ComBorders
        ComBorderNodes com_border_nodes; 	    // List of nodes on the communication borders
        Util::Vector<long> crosspoints;         // List of cross points between processes
        int rank = 0;                           // Processor rank
        MPI_Comm comm = MPI_COMM_WORLD;         // MPI Communicator
        VectorConverter vector_converter;       // L1-L2 vector converter

    public:
    	/*
    	 * Constructor for unrefined meshes
    	 *
         * m: Number of processes in vertical direction
         * n: Number of processes in horizontal direction
         * comm: MPI communicator
         * rank: MPI process rank
         *
    	 */
        Skeleton(long m, long n, MPI_Comm comm, int rank) :
                 com_borders(2 * n * m - n - m),
                 com_border_nodes(2 * n * m - n - m),
                 n_borders(2*n*m-n-m),
                 crosspoints((m + 1) * (n + 1)),
                 comm(comm),
                 rank(rank) {}

    	/*
    	 * Constructor for refined meshes
    	 *
         * m: Number of processes in vertical direction
         * n: Number of processes in horizontal direction
         * refine_factor: number of performed refinements
         * comm: MPI communicator
         * rank: MPI process rank
         *
    	 */
        Skeleton(long m, long n, long refine_factor, MPI_Comm comm, int rank) :
                 com_borders(2 * n * m - n - m),
                 com_border_nodes(2 * n * m - n - m, pow(2, refine_factor) - 1),
                 n_borders(2*n*m-n-m),
                 crosspoints((m + 1) * (n + 1)),
                 comm(comm), rank(rank) {}

        /*
    	 * Constructor for local meshes
    	 *
         * n_borders: Number of communication borders
         * nodes_per_border: Number of nodes on each communication border
         * comm: MPI communicator
         * rank: MPI process rank
         *
    	 */
        Skeleton(long n_borders, long nodes_per_border, MPI_Comm comm,
                 int rank, enum global_or_local usecase) :
                 com_borders(n_borders),
                 com_border_nodes(n_borders, nodes_per_border),
                 comm(comm), rank(rank),
                 n_borders(n_borders) {
#ifndef NDEBUG
            assert(usecase == LOCAL);
#endif
        }

        Skeleton(Skeleton &&) = delete;
        Skeleton(const Skeleton &) = delete;
        
        Skeleton & operator=(Skeleton &&other) = default;
        Skeleton & operator=(const Skeleton &) = delete;

        // General getter methods
        const long get_n_borders() const { return n_borders; }
        const long get_n_border_nodes() const { return com_border_nodes.get_n_nodes(); }

        const ComBorder &get_border(long ix_border) const {
        	return com_borders(ix_border);
        }

        const long get_border_node(long ix_border, long ix_node) const {
            return com_border_nodes.get_entry(ix_border, ix_node);
        }

        const Util::Vector<long> &get_crosspoints() const { return crosspoints; }

        const MPI_Comm get_comm() const { return comm; }
        const int get_rank() const { return rank; }

        void set_vector_converter(VectorConverter &&converter) { vector_converter = std::move(converter); }

    	// Creating global Skeleton from Mesh
        void Create(Mesh::Mesh &mesh);
        
    	// Transforms global Skeleton to local Skeleton
        void CreateLocal(Mesh::LocalMesh &local_mesh);

        // Scatter Skeleton between Processes by MPI
        void Scatter(Mesh::LocalMesh &local_mesh);

        // Vector transformations
        void AccumulatedToDistributed(Util::Vector<double> &local_vector) const {
            vector_converter.AccumulatedToDistributed(local_vector);
        }
        void DistributedToAccumulated(Util::Vector<double> &local_vector) const {
            vector_converter.DistributedToAccumulated(local_vector, *this);
        }
        void DistributedToAccumulated(Util::Vector<double> &local_vector_send,
                                      Util::Vector<double> &local_vector_recv) const {
            vector_converter.DistributedToAccumulated(local_vector_send, local_vector_recv, *this);
        }
        void GatherAccumulatedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv) const {
            vector_converter.GatherAccumulatedVector(local_vector_send, global_vector_recv, *this);
        }
        void GatherDistributedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv) const {
            vector_converter.GatherDistributedVector(local_vector_send, global_vector_recv, *this);
        }

		// Print data of Skeleton
        void Print();
    };
}

#endif //_MPI

#endif //HPC2_SKELETON_HPP
