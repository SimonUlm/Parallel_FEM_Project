#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <algorithm>
#include <array>
#include <cassert>
#ifdef _MPI
#include <mpi.h>
#endif

#include "hpc.hpp"

namespace Mesh{
    /**
     * @brief Abstract Mesh Class
     */
    class Mesh {
    protected:
		/**
	 	* @brief Number of rectangles for this mesh
		*/
		long m, n;
		long refine_factor = 0;

		Util::Vector<Node> nodes;
        Util::Vector<Element> elements;
        Util::Vector<Edge> edges;
        Util::Vector<BoundaryEdge> boundary;
        Util::Vector<long> fixed_nodes;

        // Defined in Collect.cpp
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
        
        Mesh(long m, long n, long nnodes, long nelem, long nbdry, long nedges, long nfixed) :
            m(m), n(n),
            nodes(nnodes), elements(nelem), edges(nedges),
            boundary(nbdry), fixed_nodes(nfixed) {}

        
        Mesh(Mesh &&) = delete;
        Mesh(const Mesh &) = delete;

        Mesh & operator=(Mesh &&other) = default;
        Mesh & operator=(const Mesh &) = delete;

        // Defined in Print.cpp
        void Print();

        // Defined in Collect.cpp
        void CollectEdges();

        long get_m() {return m;}
        long get_n() {return n;}
        long get_refine_factor() {return refine_factor;}

    };

    class LocalMesh: public Mesh {
    private:
        Util::VectorConverter vector_converter_;

    public:
        friend class GlobalMesh;

        Util::Vector<long> local_to_global;

        long get_n_nodes() {
            return nodes.count();
        };

        const Util::VectorConverter & vector_converter() {
            return vector_converter_;
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

        // Defined in Refine.cpp
        void Refine();

        void CollectFixedNodes() {
            Mesh::CollectFixedNodes(boundary.count());
        }

        // Defined in Scatter.cpp
#ifdef _MPI
        void Scatter(LocalMesh &local_mesh, MPI_Comm comm, int rank, int nof_processes);
    private:
        void TransferGlobalToLocal(LocalMesh &local_mesh, Util::Vector<long> &global_nodes_priority,
                                   MPI_Comm comm, int rank);
#endif
    };
}

#endif //HPC2_MESH_HPP
