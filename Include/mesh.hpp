#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <algorithm>
#include <array>

#ifdef _MPI

#include <mpi.h>

namespace Skeleton {
    class Skeleton;
}
#endif //_MPI

#include "hpc.hpp"
#include "mesh_objects.hpp"


namespace Mesh {

    /*
     * Abstract Parent Mesh class
     *
     * The abstract mesh class implements general functionality which is needed
     * by local and global meshes. It handles all data structures needed by the
     * children classes as Util::Vectors and implements the general constructors
     *
     */
    class Mesh {
    public:
        Mesh() :
                m_(0), n_(0) {}

        explicit Mesh(std::array<long, 7> mesh_data) :
                Mesh(mesh_data[0], mesh_data[1],
                     mesh_data[2], mesh_data[3], mesh_data[4],
                     mesh_data[5], mesh_data[6]) {}

        Mesh(long m, long n, long nnodes, long nelem, long nbdry) :
                Mesh(m, n, nnodes, nelem, nbdry, 0, 0) {}

        Mesh(long m, long n, long nnodes, long nelem, long nbdry,
             long nedges, long nfixed) :
                m_(m), n_(n),
                nodes_(nnodes), elements_(nelem), edges_(nedges),
                boundary_(nbdry), fixed_nodes_(nfixed) {}


        Mesh(Mesh &&) = delete;

        Mesh(const Mesh &) = delete;

        Mesh &operator=(Mesh &&) = default;

        Mesh &operator=(const Mesh &) = delete;

        // General getter methods
        const long m() const { return m_; }

        const long n() const { return n_; }

        const long n_nodes() const { return nodes_.count(); }

        const long refine_factor() const { return refine_factor_; }

        // Print data of mesh
        void Print();

        // Automatically generate a list of edges after a refinement
        void CollectEdges();

        // Generate stiffness matrix from mesh
        Util::SedMatrix CreateStiffness();

        // Generate rhs of equation system
        Util::BlasVector CreateRhs(double (*fvol)(Node &, long),
                                   double (*fNeu)(Node &, long));

        // Adding dirichlet boundary conditions to equation system
        void AddDirichlet(Util::SedMatrix &stiff_matrix,
                          Util::BlasVector &b,
                          double (*fDir)(Node &, long));

    protected:
        long m_, n_;                // Number of processors in X/Y-Direction
        long refine_factor_ = 0;    // Number of refinements performed

        Util::Vector<Node> nodes_;               // Vector of Nodes that make up the mesh
        Util::Vector<Element> elements_;         // Vector of Elements that make up the mesh
        Util::Vector<Edge> edges_;               // Vector of Edges that make up the mesh
        Util::Vector<BoundaryEdge> boundary_;    // Vector of Edges on the Boundary of the mesh
        Util::Vector<long> fixed_nodes_;         // Vector of Nodes that are fixed

        // Automatically generate a list of fixed nodes_ after a refinement
        void CollectFixedNodes(long global_nbdry);
    };


    /*
     * LocalMesh
     *
     * Implements the special functionality for local meshes. It saves the mapping
     * between the global and local indices in a vector
     *
     */
    class LocalMesh : public Mesh {
    public:
        // Befriend LocalMesh class so each other can access private members
        friend class GlobalMesh;

        // Mapping between local and global indices
        Util::Vector<long> local_to_global;

        /*
         * Constructor inherited from parent class
         *
         */
        LocalMesh() :
                Mesh::Mesh() {}

        /*
         * Constructor inherited from parent class.
         * Also initializes local_to_global
         *
         */
        explicit LocalMesh(std::array<long, 7> mesh_data) :
                Mesh::Mesh(mesh_data),
                local_to_global(mesh_data[2]) {}

        LocalMesh(LocalMesh &&) = delete;

        LocalMesh(const LocalMesh &) = delete;

        LocalMesh &operator=(LocalMesh &&) = default;

        LocalMesh &operator=(const LocalMesh &) = delete;

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(0);
        }
    };


    /*
     * GlobalMesh
     *
     * Implements the special functionality for global meshes. Such as the initial
     * creation, the refinement and the scattering to the local processes.
     *
     */
    class GlobalMesh : public Mesh {
    public:
        // Befriend LocalMesh class so each other can access private members
        friend class LocalMesh;

        /*
         * Constructor inherited from parent class
         *
         */
        GlobalMesh() :
                Mesh::Mesh() {}

        /*
         * Constructor inherited from parent class
         *
         */
        explicit GlobalMesh(std::array<long, 7> mesh_data) :
                Mesh::Mesh(mesh_data) {}

        /*
         * Constructor inherited from parent class
         *
         */
        GlobalMesh(long m, long n, long nnodes, long nelem, long nbdry) :
                Mesh::Mesh(m, n, nnodes, nelem, nbdry) {}

        /*
         * General constructor which calculates number of elements_ from number of
         * processes
         *
         */
        GlobalMesh(long m, long n) :
                Mesh(m, n,
                     (m + 1) * (n + 1), // nnodes
                     2 * m * n,         // nelem
                     2 * (m + n),       // nbdrdy
                     3 * m * n + m + n, // nedges
                     m + n + 1) {}      // nfixed

        GlobalMesh(GlobalMesh &&) = delete;

        GlobalMesh(const GlobalMesh &) = delete;

        GlobalMesh &operator=(GlobalMesh &&) = default;

        GlobalMesh &operator=(const GlobalMesh &) = delete;

        // Create rectangular unrefined global mesh
        void Create(Node bottom_left_node = Node{0, 0}, Node top_right_node = Node{1, 1});

        void Create(long m, long n, Node bottom_left_node = Node{0, 0}, Node top_right_node = Node{1, 1}) {
            *this = GlobalMesh(m, n);
            Create(bottom_left_node, top_right_node);
        }

        // Refine global mesh
        void Refine();

        void Refine(int refine_factor) {
            for (int i = 0; i < refine_factor; ++i)
                Refine();
        }

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(boundary_.count());
        }

        // Defined in mesh_scatter.cpp
#ifdef _MPI

        void Scatter(LocalMesh &local_mesh, Skeleton::Skeleton &skeleton);

    private:
        void TransferGlobalToLocal(LocalMesh &local_mesh,
                                   Util::Vector<long> &global_nodes_priority,
                                   MPI_Comm comm, int rank);

#endif
    };
}

#endif //HPC2_MESH_HPP