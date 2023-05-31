#include <memory>

#include "hpc.hpp"

namespace Mesh {

    void Mesh::CollectEdges() {

        // Get the number of edges
        long nof_edges = 0;
        for (auto element : elements)
            for (long k = 1; k <= 3; ++k)
                if (element.get_m(k) > nof_edges)
                    nof_edges = element.get_m(k);

        // Allocate storage for edge information (+ 1 because we started counting at 0)
        edges = Util::Vector<Edge>(nof_edges + 1);

        // Get endpoints for each edge, i.e. compute edgeno for the edges
        for (auto element : elements)
            for (long k = 1; k <= 3; ++k)
                edges(element.get_m(k)) = Edge{
                        element.get_n(k), // n1
                        element.get_successor_n(k) // n2
                };
    }

    void Mesh::CollectFixedNodes(long global_nbdry) {

        // Use the affiliation to set a flag for each dirichlet node
        std::unique_ptr<bool[]> flags(new bool[nodes.count()]{});
        long count_fixed = 0;
        for (auto boundary_edge : boundary) {
            if (boundary_edge.t == 1) {
                flags[boundary_edge.n1] = true;
                flags[boundary_edge.n2] = true;
                ++count_fixed;
            }
        }

        // Number of fixed nodes is higher by one than number of fixed boundary edges
        // (Except the mesh is global and the whole boundary is fixed)
        if (count_fixed != 0 && count_fixed != global_nbdry)
            ++count_fixed;

        // Write to fixed nodes
        fixed_nodes = Util::Vector<long>(count_fixed);
        long index = 0;
        for (long j = 0; j < nodes.count(); ++j)
            if (flags[j])
                fixed_nodes(index++) = j;
    }
}
