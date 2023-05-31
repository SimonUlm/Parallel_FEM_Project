#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

namespace Mesh {

    void GlobalMesh::Scatter(LocalMesh &local_mesh, MPI_Comm comm, int rank, int nof_processes) {

        // Declare temporary global mesh structure that is used by all processes. This approach is chosen to make sure
        // the global mesh gets destructed at the end of the scatter method on all local processes to save memory.
        GlobalMesh global_mesh_temp;

        // Allocate global mesh on each process
        std::array<long, 7> mesh_data{};
        if (rank == 0)
            mesh_data = {m, n, nodes.count(), elements.count(), boundary.count(), edges.count(), fixed_nodes.count()};
        MPI_Bcast(mesh_data.data(), 7, MPI_LONG, 0, comm);
        if (rank == 0)
            global_mesh_temp = std::move(*this);
        else
            global_mesh_temp = GlobalMesh(mesh_data);

        // Send global mesh to all processes
        MPI_Bcast(global_mesh_temp.nodes.data(), (int) global_mesh_temp.nodes.count() * 2,
                  MPI_DOUBLE, 0, comm);
        MPI_Bcast(global_mesh_temp.elements.data(), (int) global_mesh_temp.elements.count() * 7,
                  MPI_LONG, 0, comm);
        MPI_Bcast(global_mesh_temp.boundary.data(), (int) global_mesh_temp.boundary.count() * 4,
                  MPI_LONG, 0, comm);

        // Prepare data structure that counts for all nodes by how many processes shared with
        Util::Vector<long> global_nodes_priority(global_mesh_temp.nodes.count());

        // Gather relevant information and write into local mesh
        global_mesh_temp.TransferGlobalToLocal(local_mesh, global_nodes_priority, comm, rank);
        local_mesh.CollectEdges();
        local_mesh.CollectFixedNodes();

        // Create vector converter
        local_mesh.vector_converter_ = Util::VectorConverter(global_nodes_priority,
                                                             local_mesh.local_to_global);

        // Make sure the root process gets its global mesh back
        if (rank == 0)
            *this = std::move(global_mesh_temp);
    }

    void GlobalMesh::TransferGlobalToLocal(LocalMesh &local_mesh, Util::Vector<long> &global_nodes_priority,
                                           MPI_Comm comm, int rank) {

        // Create four temporary arrays that maps global to local nodes
        // 1. Flag array that determines whether a global edge is also local
        std::unique_ptr<long[]> edge_flags(new long[edges.count()]{});
        // 2. Flag array that determines whether a global node is also local
        std::unique_ptr<long[]> node_flags(new long[nodes.count()]{});
        // 3. Array that maps global to local edges
        std::unique_ptr<long[]> global_to_local_edges(new long[edges.count()]{});
        // 4. Array that maps global to local nodes
        std::unique_ptr<long[]> global_to_local_nodes(new long[nodes.count()]{});

        // Count elements and remember the global nodes and boundary edges that belong to the local mesh
        int elem_count = 0;
        for (auto element : elements) {
            if (element.t == rank) {
                node_flags[element.n1] = true;
                node_flags[element.n2] = true;
                node_flags[element.n3] = true;
                edge_flags[element.m1] = true;
                edge_flags[element.m2] = true;
                edge_flags[element.m3] = true;
                ++elem_count;
            }
        }

        // Create temporary global_to_local_edges
        int edge_count = 0;
        for (long i = 0; i < edges.count(); ++i) {
            if (edge_flags[i]) {
                global_to_local_edges[i] = edge_count;
                ++edge_count;
            }
        }

        // Create temporary global_to_local_nodes
        int node_count = 0;
        for (long i = 0; i < nodes.count(); ++i) {
            if (node_flags[i]) {
                global_to_local_nodes[i] = node_count;
                ++node_count;
            }
        }

        // Create local_to_global member
        local_mesh.local_to_global = Util::Vector<long>(node_count);
        node_count = 0;
        for (long i = 0; i < nodes.count(); ++i) {
            if (node_flags[i]) {
                local_mesh.local_to_global(node_count) = i;
                ++node_count;
            }
        }

        // Transfer and renumber elements
        local_mesh.elements = Util::Vector<Element>(elem_count);
        elem_count = 0;
        for (auto element : elements) {
            if (element.t == rank) {
                local_mesh.elements(elem_count) = element;
                local_mesh.elements(elem_count).n1 = global_to_local_nodes[element.n1];
                local_mesh.elements(elem_count).n2 = global_to_local_nodes[element.n2];
                local_mesh.elements(elem_count).n3 = global_to_local_nodes[element.n3];
                local_mesh.elements(elem_count).m1 = global_to_local_edges[element.m1];
                local_mesh.elements(elem_count).m2 = global_to_local_edges[element.m2];
                local_mesh.elements(elem_count).m3 = global_to_local_edges[element.m3];
                ++elem_count;
            }
        }

        // Count boundary edges
        int bdry_count = 0;
        for (auto boundary_edge : boundary)
            if (edge_flags[boundary_edge.m])
                ++bdry_count;

        // Transfer and renumber boundary edges
        if (bdry_count != 0) {
            local_mesh.boundary = Util::Vector<BoundaryEdge>(bdry_count);
            bdry_count = 0;
            for (auto boundary_edge : boundary) {
                if (edge_flags[boundary_edge.m]) {
                    local_mesh.boundary(bdry_count) = boundary_edge;
                    local_mesh.boundary(bdry_count).n1 = global_to_local_nodes[boundary_edge.n1];
                    local_mesh.boundary(bdry_count).n2 = global_to_local_nodes[boundary_edge.n2];
                    local_mesh.boundary(bdry_count).m = global_to_local_edges[boundary_edge.m];
                    ++bdry_count;
                }
            }
        }

        // Transfer nodes
        local_mesh.nodes = Util::Vector<Node>(node_count);
        for (long i = 0; i < node_count; ++i)
            local_mesh.nodes(i) = nodes(local_mesh.local_to_global(i));

        // Calculate priority for each node
        MPI_Allreduce(node_flags.get(), global_nodes_priority.data(), (int) nodes.count(),
                      MPI_LONG, MPI_SUM, comm);
    }
}
#endif // _MPI
