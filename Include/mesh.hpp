#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <algorithm>
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

        List<Node> nodes;
        List<Element> elements;
        List<Edge> edges;
        List<BoundaryEdge> boundary;
        List<long> fixed_nodes;

    public:
        Mesh() :
            m(0), n(0) {}

        explicit Mesh(std::array<long, 7> mesh_data) :
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
        void CollectFixedNodes();
        
        long get_m() {return m;}
        long get_n() {return n;}
    };

    class LocalMesh: public Mesh {
    private:
        /**
     * @brief Local Number of rectangles
    */
        long k, l;
        /**
	 * @brief Rank of Local process
	*/
        long rank;
    public:
        friend class GlobalMesh;

        List<long> local_to_global;
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

        // Defined in Scatter.cpp
#ifdef _MPI
        void Scatter(LocalMesh &local_mesh, MPI_Comm comm, int rank, int nof_processes);
    private:
        void TransferGlobalToLocal(LocalMesh &local_mesh, int rank);
#endif
    };
}



// Other declarations
#ifdef _MPI
void MpiPrintSerial(Mesh::LocalMesh &mesh, MPI_Comm comm, int rank, int nof_processes);
#endif

#endif //HPC2_MESH_HPP
