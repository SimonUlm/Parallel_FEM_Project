#include "hpc.hpp"

namespace Mesh {

    /*
     * Refine global mesh by splitting each edge into two edges_ with a new node in
     * the middle. Nodes with the same coordinates after the refinement have still
     * the same index as before. All other elements_ are completely new numbered.
     * The mesh is overwritten with the new refined mesh
     *
     */
    void GlobalMesh::Refine() {
        // Allocate new mesh
        long nnodes = nodes_.count() + edges_.count();
        long nelem = elements_.count() * 4;
        long nbdry = boundary_.count() * 2;

        GlobalMesh new_mesh(m_, n_, nnodes, nelem, nbdry);
        new_mesh.refine_factor_ = refine_factor_ + 1;

        // Copy old coordinates
        for (long i = 0; i < nodes_.count(); ++i)
            new_mesh.nodes_(i) = nodes_(i);

        // Compute new coordinates
        for (long i = 0; i < edges_.count(); ++i)
            new_mesh.nodes_(nodes_.count() + i) = (nodes_(edges_(i).n1)
                                                   + nodes_(edges_(i).n2)) * 0.5;

        // Create new elements
        for (long i = 0; i < elements_.count(); ++i) {
            // old element
            Element element = elements_(i);

            // new outer triangles
            for (long k = 1; k < 4; ++k)
                new_mesh.elements_(4 * i + k - 1) = Element{
                        element.get_n(k), // n1
                        element.get_m(k) + nodes_.count(), // n2
                        element.get_predecessor_m(k) + nodes_.count(), // n3
                        element.get_m(k) * 2 + (element.get_n(k) > element.get_successor_n(k)), // m1
                        edges_.count() * 2 + i * 3 + k - 1, // m2
                        element.get_predecessor_m(k) * 2 + (element.get_n(k) > element.get_predecessor_n(k)), // m3
                        element.t // t
                };

            // new inner triangle
            new_mesh.elements_(4 * i + 3) = Element{
                    element.m1 + nodes_.count(), // n1
                    element.m2 + nodes_.count(), // n2
                    element.m3 + nodes_.count(), // n3
                    edges_.count() * 2 + i * 3 + 1, // m1
                    edges_.count() * 2 + i * 3 + 2, // m2
                    edges_.count() * 2 + i * 3 + 0, // m3
                    element.t
            };
        }

        // Create new boundary
        for (long i = 0; i < boundary_.count(); ++i) {
            // old boundary_ edge
            BoundaryEdge edge = boundary_(i);

            // first new boundary edge
            new_mesh.boundary_(2 * i) = BoundaryEdge{
                    edge.n1, //n1
                    edge.m + nodes_.count(), // n2
                    edge.m * 2 + (edge.n1 > edge.n2), // m_
                    edge.t // t
            };

            // second new boundary edge
            new_mesh.boundary_(2 * i + 1) = BoundaryEdge{
                    edge.m + nodes_.count(), // n1
                    edge.n2, // n2
                    edge.m * 2 + (edge.n1 < edge.n2), // m_
                    edge.t // t
            };
        }

        // Automatically create list of edges and fixed nodes
        new_mesh.CollectEdges();
        new_mesh.CollectFixedNodes();

        // Write refined mesh to current object, while the coarse mesh gets destructed
        *this = std::move(new_mesh);
    }
}
