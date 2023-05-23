#include <cstdio>
#include <memory>
#include "hpc.hpp"

#ifdef _MPI
#include <mpi.h>

using namespace Mesh;

void GlobalMesh::Scatter(LocalMesh &local_mesh, MPI_Comm comm, int rank, int nof_processes) {

    // Allocate global mesh on local level
    long mesh_data[7];
    if (rank == 0) {
        mesh_data[0] = m;
        mesh_data[1] = n;
        mesh_data[2] = nodes.count;
        mesh_data[3] = elements.count;
        mesh_data[4] = boundary.count;
        mesh_data[5] = edges.count;
        mesh_data[6] = fixed_nodes.count;
    }
    MPI_Bcast(mesh_data, 7, MPI_LONG, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        m = mesh_data[0];
        n = mesh_data[1];
        nodes = List<Node>(mesh_data[2]);
        elements = List<Element>(mesh_data[3]);
        boundary = List<BoundaryEdge>(mesh_data[4]);
        edges = List<Edge>(mesh_data[5]);
        fixed_nodes = List<long>(mesh_data[6]);
    }

    // Send global mesh to all nodes
    MPI_Bcast(&nodes(0).x, (int) nodes.count * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&elements(0).n1, (int) elements.count * 7, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boundary(0).n1, (int) boundary.count * 4, MPI_LONG, 0, MPI_COMM_WORLD);

    // Gather relevant information and write into local mesh
    CollectLocalElements(local_mesh, rank);
}

void GlobalMesh::CollectLocalElements(LocalMesh &local_mesh, int rank) {

    // Create three temporary arrays that maps global to local nodes
    // 1. Flag array that determines whether a global edge is also local
    std::unique_ptr<bool[]> edge_flags(new bool[edges.count]{});
    // 2. Flag array that determines whether a global node is also local
    std::unique_ptr<bool[]> node_flags(new bool[nodes.count]{});
    // 3. Array that maps global to local nodes
    std::unique_ptr<long[]> global_to_local(new long[nodes.count]{});

    // Count elements and remember the global nodes and boundary edges that belong to the local mesh
    int elem_count = 0;
    for (long i = 0; i < elements.count; ++i) {
        if (elements(i).t == rank) {
            edge_flags[elements(i).m1] = true;
            edge_flags[elements(i).m2] = true;
            edge_flags[elements(i).m3] = true;
            node_flags[elements(i).n1] = true;
            node_flags[elements(i).n2] = true;
            node_flags[elements(i).n3] = true;
            ++elem_count;
        }
    }

    // Create temporary global-to-local array
    int node_count = 0;
    for (long i = 0; i < nodes.count; ++i) {
        if (node_flags[i]) {
            global_to_local[i] = node_count;
            ++node_count;
        }
    }

    // Create local-to-global array
    local_mesh.local_to_global = List<long>(node_count);
    node_count = 0;
    for (long i = 0; i < nodes.count; ++i) {
        if (node_flags[i]) {
            local_mesh.local_to_global(node_count) = i;
            ++node_count;
        }
    }

    // Transfer and renumber elements
    local_mesh.elements = List<Element>(elem_count);
    elem_count = 0;
    for (long i = 0; i < elements.count; ++i) {
        if (elements(i).t == rank) {
            local_mesh.elements(elem_count) = elements(i);
            local_mesh.elements(elem_count).n1 = global_to_local[elements(i).n1];
            local_mesh.elements(elem_count).n2 = global_to_local[elements(i).n2];
            local_mesh.elements(elem_count).n3 = global_to_local[elements(i).n3];
            ++elem_count;
        }
    }

    // Count boundary edges
    int bdry_count = 0;
    for (long i = 0; i < boundary.count; ++i)
        if (edge_flags[boundary(i).m])
            ++bdry_count;

    // Transfer and renumber boundary edges
    if (bdry_count != 0) {
        local_mesh.boundary = List<BoundaryEdge>(bdry_count);
        bdry_count = 0;
        for (long i = 0; i < boundary.count; ++i) {
            if (edge_flags[boundary(i).m]) {
                local_mesh.boundary(bdry_count) = boundary(i);
                local_mesh.boundary(bdry_count).n1 = global_to_local[boundary(i).n1];
                local_mesh.boundary(bdry_count).n2 = global_to_local[boundary(i).n2];
                ++bdry_count;
            }
        }
    }

    // Transfer nodes
    local_mesh.nodes = List<Node>(node_count);
    for (long i = 0; i < node_count; ++i)
        local_mesh.nodes(i) = nodes(local_mesh.local_to_global(i));
}

#endif // _MPI
