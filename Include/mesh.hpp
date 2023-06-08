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



namespace Mesh{

    /*
     * Abstract Parent Mesh class
     *
     * The abstract mesh class implements general functionality which is needed
     * by local and global meshes. It handles all data structures needed by the
     * children classes as Util::Vector's and implements the general constructors
     *
     */
    class Mesh {
    protected:
		long m, n;                 // Number of processors in X/Y-Direction
		long refine_factor = 0;    // Number of refinements performed

		Util::Vector<Node> nodes;               // Vector of Nodes that make up the mesh
        Util::Vector<Element> elements;         // Vector of Elements that make up the mesh
        Util::Vector<Edge> edges;               // Vector of Edges that make up the mesh
        Util::Vector<BoundaryEdge> boundary;    // Vector of Edges on the Boundary of the mesh
        Util::Vector<long> fixed_nodes;         // Vector of Nodes that are fixed

        // Automatically generate a list of fixed nodes after a refinement
        void CollectFixedNodes(long global_nbdry);

    public:    
        Mesh() :
                m(0), n(0) {}

        Mesh(std::array<long, 7> mesh_data) :
                Mesh(mesh_data[0], mesh_data[1],
                 mesh_data[2], mesh_data[3], mesh_data[4],
                 mesh_data[5], mesh_data[6]) {}

        Mesh(long m, long n, long nnodes, long nelem, long nbdry) :
                Mesh(m, n, nnodes, nelem, nbdry, 0, 0) {}
        
        Mesh(long m, long n, long nnodes, long nelem, long nbdry,
             long nedges, long nfixed) :
             m(m), n(n),
             nodes(nnodes), elements(nelem), edges(nedges),
             boundary(nbdry), fixed_nodes(nfixed) {}

        
        Mesh(Mesh &&) = delete;
        Mesh(const Mesh &) = delete;

        Mesh & operator=(Mesh &&other) = default;
        Mesh & operator=(const Mesh &) = delete;

        // General getter methods
        const long get_m() const { return m; }
        const long get_n() const { return n; }
        const long get_n_nodes() const { return nodes.count(); }
        const long get_refine_factor() const { return refine_factor; }

        // Print data of mesh
        void Print();

        // Automatically generate a list of edges after a refinement
        void CollectEdges();
        
        // Generate stiffness matrix from mesh
        Util::SedMatrix CreateStiffness();

        // Generate rhs of equation system
	    Util::BlasVector CreateRhs(double (*fvol)(Node&, long),
                                   double (*fNeu)(Node&, long));

		// Adding dirichlet boundary conditions to equation system
	    void AddDirichlet(Util::SedMatrix &stiff_matrix,
			              Util::BlasVector &b,
			              double (*fDir)(Node&, long));
    };


    /*
     * LocalMesh
     *
     * Implements the special functionality for local meshes. It saves the mapping
     * between the global and local indices in a vector
     *
     */
    class LocalMesh: public Mesh {
    public:
        // Befriend LocalMesh class so each other can access private members
        friend class GlobalMesh;

        // Mapping between local and global indices
        Util::Vector<long> local_to_global;

        LocalMesh():
                Mesh::Mesh() {}

        LocalMesh(std::array<long, 7> mesh_data) :
                Mesh::Mesh(mesh_data),
                local_to_global(mesh_data[2]) {}

        long get_n_nodes() {
            return nodes.count();
        };

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(0);
        }
    };


    /*
     * GlobalMesh
     *
     * Implements the special functionality for global meshes. Such as the inital
     * creation, the refinement and the scattering to the local processes.
     *
     */
    class GlobalMesh: public Mesh {
    public:
        using Mesh::Mesh;

        // Befriend LocalMesh class so each other can access private members
        friend class LocalMesh;

        /*
         * General constructor which calculates number of elements from number of
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
            Mesh::CollectFixedNodes(boundary.count());
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
