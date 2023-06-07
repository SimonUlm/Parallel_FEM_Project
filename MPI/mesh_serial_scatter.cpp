#include <cstdio>
#include <memory>
#include <array>
#include "hpc.hpp"

#ifdef _MPI

#include <mpi.h>

#ifdef _SERIAL_SCATTER
namespace Mesh {

    void GlobalMesh::Scatter(LocalMesh &local_mesh, Skeleton::Skeleton &skeleton) {

        // MPI stuff
        MPI_Comm comm = skeleton.get_comm();
        int rank = skeleton.get_rank();
        long nnodes;
        if (rank == 0)
            nnodes = nodes.count();
        MPI_Bcast(&nnodes, 1, MPI_LONG, 0, comm);

        // Temporary
        LocalMesh mesh_sender;
        Util::Vector<long> global_nodes_priority(nnodes);
        std::array<long, 7> mesh_data{};
        MPI_Request request;
        MPI_Status status;

        if (rank == 0) {
            for (int process = 1; process < m * n; ++process) {
                local_mesh.m = m;
                local_mesh.n = n;
                TransferGlobalToLocal(local_mesh, global_nodes_priority, comm, process);
                if (process != 1)
                    MPI_Wait(&request, &status);
                mesh_sender = std::move(local_mesh);
                mesh_data = {mesh_sender.m, mesh_sender.n,
                             mesh_sender.nodes.count(), mesh_sender.elements.count(), mesh_sender.boundary.count(),
                             mesh_sender.edges.count(), mesh_sender.fixed_nodes.count()};
                MPI_Isend(mesh_data.data(), 7,
                          MPI_LONG, process, 0, comm, &request);
                MPI_Isend(mesh_sender.nodes.data(), (int) mesh_sender.nodes.count() * 2,
                         MPI_DOUBLE, process, 0, comm, &request);
                MPI_Isend(mesh_sender.elements.data(), (int) mesh_sender.elements.count() * 7,
                         MPI_LONG, process, 0, comm, &request);
                MPI_Isend(mesh_sender.boundary.data(), (int) mesh_sender.boundary.count() * 4,
                         MPI_LONG, process, 0, comm, &request);
                MPI_Isend(mesh_sender.local_to_global.data(), (int) mesh_sender.local_to_global.count(),
                          MPI_LONG, process, 0, comm, &request);
            }
            TransferGlobalToLocal(local_mesh, global_nodes_priority, comm, rank);
        } else {
            MPI_Recv(mesh_data.data(), 7, MPI_LONG, 0, 0, comm, &status);
            local_mesh = LocalMesh(mesh_data);
            MPI_Recv(local_mesh.nodes.data(), (int) local_mesh.nodes.count() * 2,
                      MPI_DOUBLE, 0, 0, comm, &status);
            MPI_Recv(local_mesh.elements.data(), (int) local_mesh.elements.count() * 7,
                      MPI_LONG, 0, 0, comm, &status);
            MPI_Recv(local_mesh.boundary.data(), (int) local_mesh.boundary.count() * 4,
                      MPI_LONG, 0, 0, comm, &status);
            MPI_Recv(local_mesh.local_to_global.data(), (int) local_mesh.local_to_global.count(),
                      MPI_LONG, 0, 0, comm, &status);
        }
        MPI_Bcast(global_nodes_priority.data(), (int) global_nodes_priority.count(), MPI_LONG, 0, comm);

        local_mesh.CollectEdges();
        local_mesh.CollectFixedNodes();

        // Scatter skeleton
        skeleton.Scatter(local_mesh);

        // Create vector converter
        skeleton.set_vector_converter(Skeleton::VectorConverter(mesh_data[0], mesh_data[1],
                                                                global_nodes_priority, local_mesh.local_to_global));
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
        for (auto &element : elements) {
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
        for (auto &element : elements) {
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
        for (auto &boundary_edge : boundary)
            if (edge_flags[boundary_edge.m])
                ++bdry_count;

        // Transfer and renumber boundary edges
        if (bdry_count != 0) {
            local_mesh.boundary = Util::Vector<BoundaryEdge>(bdry_count);
            bdry_count = 0;
            for (auto &boundary_edge : boundary) {
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
        for (long i = 0; i < nodes.count(); ++i)
            global_nodes_priority(i) += node_flags[i];
    }
}
#endif // _SERIAL_SCATTER

#endif // _MPI
