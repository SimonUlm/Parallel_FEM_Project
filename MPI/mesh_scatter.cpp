#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

namespace Mesh {

    void GlobalMesh::Scatter(LocalMesh &local_mesh, Skeleton::Skeleton &skeleton, MPI_Comm comm, int rank) {

        // Declare temporary global mesh structure that is used by all processes. This approach is chosen to make sure
        // the global mesh gets destructed at the end of the scatter method on all local processes to save memory.
        GlobalMesh global_mesh_temp;

        // Allocate global mesh on each process
        std::array<long, 7> mesh_data{};
        if (rank == 0)
            mesh_data = {m_, n_, nodes_.count(), elements_.count(), boundary_.count(), edges_.count(),
                         fixed_nodes_.count()};
        MPI_Bcast(mesh_data.data(), 7, MPI_LONG, 0, comm);
        if (rank == 0)
            global_mesh_temp = std::move(*this);
        else
            global_mesh_temp = GlobalMesh(mesh_data);

        // Send global mesh to all processes
        MPI_Bcast(global_mesh_temp.nodes_.data(), (int) global_mesh_temp.nodes_.count() * 2,
                  MPI_DOUBLE, 0, comm);
        MPI_Bcast(global_mesh_temp.elements_.data(), (int) global_mesh_temp.elements_.count() * 7,
                  MPI_LONG, 0, comm);
        MPI_Bcast(global_mesh_temp.boundary_.data(), (int) global_mesh_temp.boundary_.count() * 4,
                  MPI_LONG, 0, comm);

        // Prepare data_ structure that counts for all nodes_ by how many processes shared with
        Util::Vector<long> global_nodes_priority(global_mesh_temp.nodes_.count());

        // Gather relevant information and write into local mesh
        global_mesh_temp.TransferGlobalToLocal(local_mesh, global_nodes_priority, comm, rank);
        local_mesh.CollectEdges();
        local_mesh.CollectFixedNodes();

        // Scatter skeleton
        long refine_factor = refine_factor_;
        MPI_Bcast(&refine_factor, 1, MPI_LONG, 0, comm);
        skeleton.Scatter(local_mesh, mesh_data[0], mesh_data[1], refine_factor, comm, rank);

        // Create vector converter
        skeleton.set_vector_converter(Skeleton::VectorConverter(mesh_data[0], mesh_data[1],
                                                                global_nodes_priority, local_mesh.local_to_global));

        // Make sure the root process gets its global mesh back
        if (rank == 0)
            *this = std::move(global_mesh_temp);
    }

    void GlobalMesh::TransferGlobalToLocal(LocalMesh &local_mesh, Util::Vector<long> &global_nodes_priority,
                                           MPI_Comm comm, int rank) {

        // Create four temporary arrays that maps global to local nodes_
        // 1. Flag array that determines whether a global edge is also local
        std::unique_ptr<long[]> edge_flags(new long[edges_.count()]{});
        // 2. Flag array that determines whether a global node is also local
        std::unique_ptr<long[]> node_flags(new long[nodes_.count()]{});
        // 3. Array that maps global to local edges_
        std::unique_ptr<long[]> global_to_local_edges(new long[edges_.count()]{});
        // 4. Array that maps global to local nodes_
        std::unique_ptr<long[]> global_to_local_nodes(new long[nodes_.count()]{});

        // Count elements_ and remember the global nodes_ and boundary_ edges_ that belong to the local mesh
        int elem_count = 0;
        for (auto &element: elements_) {
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
        for (long i = 0; i < edges_.count(); ++i) {
            if (edge_flags[i]) {
                global_to_local_edges[i] = edge_count;
                ++edge_count;
            }
        }

        // Create temporary global_to_local_nodes
        int node_count = 0;
        for (long i = 0; i < nodes_.count(); ++i) {
            if (node_flags[i]) {
                global_to_local_nodes[i] = node_count;
                ++node_count;
            }
        }

        // Create local_to_global member
        local_mesh.local_to_global = Util::Vector<long>(node_count);
        node_count = 0;
        for (long i = 0; i < nodes_.count(); ++i) {
            if (node_flags[i]) {
                local_mesh.local_to_global(node_count) = i;
                ++node_count;
            }
        }

        // Transfer and renumber elements_
        local_mesh.elements_ = Util::Vector<Element>(elem_count);
        elem_count = 0;
        for (auto &element: elements_) {
            if (element.t == rank) {
                local_mesh.elements_(elem_count) = element;
                local_mesh.elements_(elem_count).n1 = global_to_local_nodes[element.n1];
                local_mesh.elements_(elem_count).n2 = global_to_local_nodes[element.n2];
                local_mesh.elements_(elem_count).n3 = global_to_local_nodes[element.n3];
                local_mesh.elements_(elem_count).m1 = global_to_local_edges[element.m1];
                local_mesh.elements_(elem_count).m2 = global_to_local_edges[element.m2];
                local_mesh.elements_(elem_count).m3 = global_to_local_edges[element.m3];
                ++elem_count;
            }
        }

        // Count boundary_ edges_
        int bdry_count = 0;
        for (auto &boundary_edge: boundary_)
            if (edge_flags[boundary_edge.m])
                ++bdry_count;

        // Transfer and renumber boundary_ edges_
        if (bdry_count != 0) {
            local_mesh.boundary_ = Util::Vector<BoundaryEdge>(bdry_count);
            bdry_count = 0;
            for (auto &boundary_edge: boundary_) {
                if (edge_flags[boundary_edge.m]) {
                    local_mesh.boundary_(bdry_count) = boundary_edge;
                    local_mesh.boundary_(bdry_count).n1 = global_to_local_nodes[boundary_edge.n1];
                    local_mesh.boundary_(bdry_count).n2 = global_to_local_nodes[boundary_edge.n2];
                    local_mesh.boundary_(bdry_count).m = global_to_local_edges[boundary_edge.m];
                    ++bdry_count;
                }
            }
        }

        // Transfer nodes_
        local_mesh.nodes_ = Util::Vector<Node>(node_count);
        for (long i = 0; i < node_count; ++i)
            local_mesh.nodes_(i) = nodes_(local_mesh.local_to_global(i));

        // Calculate priority for each node
        MPI_Allreduce(node_flags.get(), global_nodes_priority.data(), (int) nodes_.count(),
                      MPI_LONG, MPI_SUM, comm);
    }
}

#endif // _MPI
