#include "hpc.hpp"

namespace Mesh{

    void Mesh::Refine() {
        // Allocate new mesh
        long nnodes = nodes.count + edges.count;
        long nelem = elements.count * 4;
        long nbdry = boundary.count * 2;
        Mesh new_mesh(m, n, nnodes, nelem, nbdry);
        new_mesh.refine_factor = refine_factor + 1;

        // Copy old coordinates
        for (long i = 0; i < nodes.count; ++i)
            new_mesh.nodes(i) = nodes(i);

        // Compute new coordinates
        for (long i = 0; i < edges.count; ++i)
            new_mesh.nodes(nodes.count + i) = (nodes(edges(i).n1)
                    + nodes(edges(i).n2)) * 0.5;

        // Create new elements
        for (long i = 0; i < elements.count; ++i) {
            // old element
            Element element = elements(i);

            // new outer triangles
            for (long k = 1; k < 4; ++k)
                new_mesh.elements(4*i+k-1) = Element{
                    element.get_n(k), // n1
                    element.get_m(k) + nodes.count, // n2
                    element.get_predecessor_m(k) + nodes.count, // n3
                    element.get_m(k) * 2 + (element.get_n(k) > element.get_successor_n(k)), // m1
                    edges.count * 2 + i * 3 + k - 1, // m2
                    element.get_predecessor_m(k) * 2 + (element.get_n(k) > element.get_predecessor_n(k)), // m3
                    element.t // t
                };

            // new inner triangle
            new_mesh.elements(4*i+3) = Element{
                element.m1 + nodes.count, // n1
                element.m2 + nodes.count, // n2
                element.m3 + nodes.count, // n3
                edges.count * 2 + i * 3 + 1, // m1
                edges.count * 2 + i * 3 + 2, // m2
                edges.count * 2 + i * 3 + 0, // m3
                element.t
            };
        }

        // Create new boundary
        for (long i = 0; i < boundary.count; ++i) {
            // old boundary edge
            BoundaryEdge edge = boundary(i);

            // first new boundary edge
            new_mesh.boundary(2*i) = BoundaryEdge{
                    edge.n1, //n1
                    edge.m + nodes.count, // n2
                    edge.m * 2 + (edge.n1 > edge.n2), // m
                    edge.t // t
            };

            // second new boundary edge
            new_mesh.boundary(2*i+1) = BoundaryEdge{
                    edge.m + nodes.count, // n1
                    edge.n2, // n2
                    edge.m * 2 + (edge.n1 < edge.n2), // m
                    edge.t // t
            };
        }

        // Write refined mesh to current object, while the coarse mesh gets destructed
        *this = std::move(new_mesh);
    }
}
