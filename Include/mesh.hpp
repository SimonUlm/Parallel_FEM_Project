#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <algorithm>
#include <array>
#include <cassert>
#ifdef _MPI
#include <mpi.h>
#endif

#include "hpc.hpp"
#include "mesh_objects.hpp"

namespace Skeleton {
    class Skeleton;
}

namespace Mesh{
    /**
     * @brief Abstract Mesh Class
     */
    class Mesh {
    protected:
		/**
	 	* @brief Number of rectangles for this mesh
		*/
		long m, n;                  // Number of processes per row n and number of Processes per column n
		long refine_factor = 0;     // Counts how often the refine method has been used on the mesh

		Util::Vector<Node> nodes;               // Vector of Nodes that make up the mesh
        Util::Vector<Element> elements;         // Vector of Elements that make up the mesh
        Util::Vector<Edge> edges;               // Vector of Edges that make up the mesh
        Util::Vector<BoundaryEdge> boundary;    // Vector of Edges on the Boundary of the mesh
        Util::Vector<long> fixed_nodes;         // Vector of Nodes that are fixed

        // Defined in Collect.cpp
        void CollectFixedNodes(long global_nbdry);      // Returns all Nodes that are Part of a Boundary Edge

    public:    
        Mesh() :                                // Initialize "Empty"
            m(0), n(0) {}

        Mesh(std::array<long, 7> mesh_data) :   
            Mesh(mesh_data[0], mesh_data[1],
                 mesh_data[2], mesh_data[3], mesh_data[4],
                 mesh_data[5], mesh_data[6]) {}

        Mesh(long m, long n, long nnodes, long nelem, long nbdry) :
            Mesh(m, n, nnodes, nelem, nbdry, 0, 0) {}
        
        Mesh(long m, long n, long nnodes, long nelem, long nbdry, long nedges, long nfixed) :
            m(m), n(n),
            nodes(nnodes), elements(nelem), edges(nedges),
            boundary(nbdry), fixed_nodes(nfixed) {}

        
        Mesh(Mesh &&) = delete;
        Mesh(const Mesh &) = delete;

        Mesh & operator=(Mesh &&other) = default;
        Mesh & operator=(const Mesh &) = delete;

        const long get_m() const { return m; }
        const long get_n() const { return n; }
        const long get_n_nodes() const { return nodes.count(); }
        const long get_refine_factor() const { return refine_factor; }

        void Print();

        void CollectEdges();
        
        Util::SedMatrix CreateStiffness();
	    Util::BlasVector CreateRhs(double (*fvol)(Node&, long),
                                   double (*fNeu)(Node&, long));
					  
	void AddDirichlet(Util::SedMatrix &stiff_matrix,
			  Util::BlasVector &b,
			  double (*fDir)(Node&, long));
    };

    class LocalMesh: public Mesh {
    public:
        friend class GlobalMesh;

        Util::Vector<long> local_to_global;

        long get_n_nodes() {
            return nodes.count();
        };

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(0);
        }
    };

    class GlobalMesh: public Mesh {
    public:
        using Mesh::Mesh;
        friend class LocalMesh;

        GlobalMesh(long m, long n) :
                Mesh(m, n,
                 (m + 1) * (n + 1), // nnodes
                 2 * m * n,         // nelem
                 2 * (m + n),       // nbdrdy
                 3 * m * n + m + n, // nedges  
                 m + n + 1) {}      // nfixed
        
        // Defined in Create.cpp
        void Create(Node bottom_left_node = Node{0, 0}, Node top_right_node = Node{1, 1});
        void Create(long m, long n, Node bottom_left_node = Node{0, 0}, Node top_right_node = Node{1, 1}) {
            *this = GlobalMesh(m, n);
            Create(bottom_left_node, top_right_node);
        }

        // Defined in Refine.cpp
        void Refine();
        void Refine(int refine_factor) {
            for (int i = 0; i < refine_factor; ++i)
                Refine();
        }

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(boundary.count());
        }

        // Defined in Scatter.cpp
#ifdef _MPI
        void Scatter(LocalMesh &local_mesh, Skeleton::Skeleton &skeleton);
    private:
        void TransferGlobalToLocal(LocalMesh &local_mesh, Util::Vector<long> &global_nodes_priority,
                                   MPI_Comm comm, int rank);
#endif
    };
}

#endif //HPC2_MESH_HPP
