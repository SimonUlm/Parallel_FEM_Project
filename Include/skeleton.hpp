#ifndef HPC2_SKELETON_HPP
#define HPC2_SKELETON_HPP

#ifdef _MPI // only needed for parallel algorithms

#include <cstdio>
#include <cmath>
#include <mpi.h>

#include "hpc.hpp"

namespace Skeleton {

    constexpr long kNumberOfColors = 4;     // number of possible colors

    enum global_or_local {
        GLOBAL, LOCAL
    }; // Enum for local skeleton creation


    /*
     * ComBorder
     *
     * The ComBorder class contains the two cross points for this border and the
     * corresponding processes. Additionally, a color is saved so each border of
     * process has a different color
     *
     */
    class ComBorder {
    public:
        // Set entries of ComBorder
        void set_entries(long i, long start_node, long end_node,
                         long left_proc, long right_proc, long couple_color) {
            index_ = i;
            c1_ = start_node;
            c2_ = end_node;
            L_ = left_proc;
            R_ = right_proc;
            color_ = couple_color;
        }

        // Copy entries from another ComBorder
        void copy_entries(ComBorder &border) const {
            border.set_entries(index_, c1_, c2_, L_, R_, color_);
        }

        // General getter methods
        const long L() const { return L_; }

        const long R() const { return R_; }

        const long c1() const { return c1_; }

        const long c2() const { return c2_; }

        const long color() const { return color_; }

        // Setter for start/end node
        void set_c1(long new_c1) { c1_ = new_c1; }

        void set_c2(long new_c2) { c2_ = new_c2; }

        void set_index(long new_ix) { index_ = new_ix; }

        // Print data_ of ComBorder
        void Print();

    private:
        long index_; // index of border
        long c1_;    // starting cross point
        long c2_;    // end cross point
        long L_;     // Left or lower process
        long R_;     // Right or upper process
        long color_; // Color of to time communication
    };


    /*
     * ComBorderNodes
     *
     * The ComBorderNodes class consists one Vector of nodes_. This vector is not
     * seperated into the single borders but is handled as one consecutive Vector.
     * The separation into each border is handled by the class itself. So from the
     * outside the Borders look seperated
     *
     */
    class ComBorderNodes {
    public:
        /*
    	 * Default Constructor
         *
    	 */
        ComBorderNodes() {}

        /*
    	 * Constructor for unrefined Mesh
         *
         * n_borders_: Number of ComBorder's in this Skeleton
         *
    	 */
        explicit ComBorderNodes(long n_borders) :
                nodes_(0),
                n_nodes_(0) {}

        /*
    	 * Constructor for refined Mesh
         *
         * n_borders_: Number of ComBorder's in this Skeleton
         * n_nodes_: Number of nodes_ per ComBorder (calculated from refinement factor)
         *
    	 */
        ComBorderNodes(long n_borders, long n_nodes) :
                nodes_(n_borders * n_nodes),
                n_nodes_(n_nodes) {}

        // General getter methods
        const long n_nodes() const { return n_nodes_; }

        const Util::Vector<long> &nodes() const { return nodes_; }

        /*
         * Get/Set node for given border
         *
         * ix_border: index_ of border
         * ix_node: index_ of node in border
         *
         */
        const long get_entry(long ix_border, long ix_node) const {
            return nodes_(ix_border * n_nodes_ + ix_node);
        }

        void set_entry(long ix_border, long ix_node, long node) {
            nodes_(ix_border * n_nodes_ + ix_node) = node;
        }

        // Print data_ ComBorderNodes seperated for each border
        void Print();

    private:
        long n_nodes_;                // Number of nodes_ per ComBorder
        Util::Vector<long> nodes_;    // Consecutive list of nodes on the ComBorders
    };


    /*
     * Skeleton
     *
     * The Skeleton is used for communication between the distributed processes.
     * The class contains a vector of ComBorders which contain the information
     * of the corresponding processes and the edge nodes. Additionally, it contains
     * an instance of ComBorderNodes which stores all nodes on the communication
     * borders.
     *
     */
    class Skeleton {
    public:
        /*
         * Default Constructor
         *
         */
        Skeleton() {}

        /*
         * Constructor for unrefined meshes
         *
         * m_: Number of processes in vertical direction
         * n_: Number of processes in horizontal direction
         * comm_: MPI communicator
         * rank_: MPI process rank
         *
         */
        Skeleton(long m, long n) :
                com_borders_(2 * n * m - n - m),
                com_border_nodes_(2 * n * m - n - m),
                n_borders_(2 * n * m - n - m),
                crosspoints_((m + 1) * (n + 1)) {}

        /*
         * Constructor for refined meshes
         *
         * m_: Number of processes in vertical direction
         * n_: Number of processes in horizontal direction
         * refine_factor_: number of performed refinements
         * comm_: MPI communicator
         * rank_: MPI process rank
         *
         */
        Skeleton(long m, long n, long refine_factor) :
                com_borders_(2 * n * m - n - m),
                com_border_nodes_(2 * n * m - n - m, pow(2, refine_factor) - 1),
                n_borders_(2 * n * m - n - m),
                crosspoints_((m + 1) * (n + 1)) {}

        /*
    	 * Constructor for local meshes
    	 *
         * n_borders_: Number of communication borders
         * nodes_per_border: Number of nodes on each communication border
         * comm_: MPI communicator
         * rank_: MPI process rank
         *
    	 */
        Skeleton(long n_borders, long nodes_per_border, MPI_Comm comm,
                 int rank, enum global_or_local usecase) :
                com_borders_(n_borders),
                com_border_nodes_(n_borders, nodes_per_border),
                comm_(comm), rank_(rank),
                n_borders_(n_borders) {
#ifndef NDEBUG
            assert(usecase == LOCAL);
#endif
        }

        Skeleton(Skeleton &&) = delete;

        Skeleton(const Skeleton &) = delete;

        Skeleton &operator=(Skeleton &&other) = default;

        Skeleton &operator=(const Skeleton &) = delete;

        // General getter methods
        const long n_borders() const { return n_borders_; }

        const long n_border_nodes() const { return com_border_nodes_.n_nodes(); }

        const Util::Vector<long> &crosspoints() const { return crosspoints_; }

        const MPI_Comm comm() const { return comm_; }

        const int rank() const { return rank_; }

        const ComBorder &get_border(long ix_border) const {
            return com_borders_(ix_border);
        }

        const long get_border_node(long ix_border, long ix_node) const {
            return com_border_nodes_.get_entry(ix_border, ix_node);
        }

        void set_vector_converter(VectorConverter &&converter) { vector_converter_ = std::move(converter); }

        // Creating global Skeleton from Mesh
        void Create(Mesh::GlobalMesh &mesh);

        void Create(Mesh::GlobalMesh &mesh, long m, long n) {
            *this = Skeleton(m, n, mesh.refine_factor());
            Create(mesh);
        }

        // Transforms global Skeleton to local Skeleton
        void CreateLocal(Mesh::LocalMesh &local_mesh, MPI_Comm comm, int rank);

        // Scatter Skeleton between Processes by MPI
        void Scatter(Mesh::LocalMesh &local_mesh, long m, long n, long refine_factor, MPI_Comm comm, int rank);

        // Vector transformations
        void AccumulatedToDistributed(Util::Vector<double> &local_vector) const {
            vector_converter_.AccumulatedToDistributed(local_vector);
        }

        void DistributedToAccumulated(Util::Vector<double> &local_vector) const {
            vector_converter_.DistributedToAccumulated(local_vector, *this);
        }

        void DistributedToAccumulated(Util::Vector<double> &local_vector_send,
                                      Util::Vector<double> &local_vector_recv) const {
            vector_converter_.DistributedToAccumulated(local_vector_send, local_vector_recv, *this);
        }

        void GatherAccumulatedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv) const {
            vector_converter_.GatherAccumulatedVector(local_vector_send, global_vector_recv, *this);
        }

        void GatherDistributedVector(Util::Vector<double> &local_vector_send,
                                     Util::Vector<double> &global_vector_recv) const {
            vector_converter_.GatherDistributedVector(local_vector_send, global_vector_recv, *this);
        }

        // Print data_ of Skeleton
        void Print();

    private:
        long n_borders_;                         // Number of borders in this Skeleton
        Util::Vector<ComBorder> com_borders_;    // Vector of ComBorders
        ComBorderNodes com_border_nodes_;        // List of nodes on the communication borders
        Util::Vector<long> crosspoints_;         // List of cross points between processes
        int rank_ = 0;                           // Processor rank
        MPI_Comm comm_ = MPI_COMM_WORLD;         // MPI Communicator
        VectorConverter vector_converter_;       // Vector converter
    };
}

#endif //_MPI

#endif //HPC2_SKELETON_HPP