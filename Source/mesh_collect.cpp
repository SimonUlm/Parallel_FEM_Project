#include <memory>
#include "hpc.hpp"

namespace Mesh {

    void Mesh::CollectEdges() {

        // Get the number of edges
        long nof_edges = 0;
        for (long i = 0; i < elements.count; ++i)
            for (long k = 1; k <= 3; ++k)
                if (elements(i).get_m(k) > nof_edges)
                    nof_edges = elements(i).get_m(k);

        // Allocate storage for edge information (+ 1 because we started counting at 0)
        edges = Util::List<Edge>(nof_edges + 1);

        // Get endpoints for each edge, i.e. compute edgeno for the edges
        for (long i = 0; i < elements.count; ++i)
            for (long k = 1; k <= 3; ++k)
                edges(elements(i).get_m(k)) = Edge{
                        elements(i).get_n(k), // n1
                        elements(i).get_successor_n(k) // n2
                };
    }

    void Mesh::CollectFixedNodes() {
        // Use the affiliation to set a flag for each dirichlet node
        std::unique_ptr<bool[]> flags(new bool[nodes.count]{});
        long count_fixed = 0;
        for (long i = 0; i < boundary.count; ++i) {
            if (boundary(i).t == 1) {
                flags[boundary(i).n1] = true;
                flags[boundary(i).n2] = true;
                ++count_fixed;
            }
        }

        // Number of fixed nodes is higher than number of fixed boundary edges (except the whole boundary is fixed)
        if (count_fixed != boundary.count)
            ++count_fixed;

        // Write to fixed nodes
        fixed_nodes = Util::List<long>(count_fixed);
        long index = 0;
        for (long j = 0; j < nodes.count; ++j)
            if (flags[j])
                fixed_nodes(index++) = j;
        fixed_nodes.count = index;
    }
}