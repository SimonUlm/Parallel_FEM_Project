#include "hpc.hpp"

namespace Mesh {

    void GlobalMesh::Create(Node bottom_left_node, Node top_right_node) {
        long node_index = 0;
        long fixed_index = 0;

        long nof_h_edges = (m + 1) * n; // number of horizontal edges
        long nof_v_edges = m * (n + 1); // number of vertical edges
        long nodes_per_row = n + 1;

        double delta_x = top_right_node.x - bottom_left_node.x;
        double delta_y = top_right_node.y - bottom_left_node.y;

        // Iterating over rectangles (2 elements at a time)
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j < n; ++j) {

                // Coordinates of bottom left node
                if (i == 0 && j == 0) {
                    // Write bottom left node of rectangle for root rectangle
                    nodes(node_index) = bottom_left_node;
                    fixed_nodes(fixed_index++) = node_index;
                }

                if (j == 0) {
                    // Write top left node for first column
                    node_index = nodes_per_row * (i + 1);
                    nodes(node_index).x = bottom_left_node.x;
                    nodes(node_index).y = bottom_left_node.y + ((double) (i + 1) / m) * delta_y;
                    fixed_nodes(fixed_index++) = node_index;
                }

                if (i == 0) {
                    // Write bottom right node for first row
                    node_index = j + 1;
                    nodes(node_index).x = bottom_left_node.x + ((double) (j + 1) / n) * delta_x;
                    nodes(node_index).y = bottom_left_node.y;
                    fixed_nodes(fixed_index++) = node_index;
                }

                // Write upper right node for each rectangle
                node_index = nodes_per_row * (i + 1) + (j + 1);
                nodes(node_index).x = bottom_left_node.x + ((double) (j + 1) / n) * delta_x;
                nodes(node_index).y = bottom_left_node.y + ((double) (i + 1) / m) * delta_y;


                // Four nodes numbers
                long node_bl1 = i * nodes_per_row + j;            // bottom left node
                long node_br2 = i * nodes_per_row + j + 1;        // bottom right node
                long node_tr3 = (i + 1) * nodes_per_row + j + 1;    // top right node
                long node_tl4 = (i + 1) * nodes_per_row + j;        // top left node

                // Five edge numbers
                long edge_bottom = i * n + j;
                long edge_top = (i + 1) * n + j;
                long edge_left = nof_h_edges + i * nodes_per_row + j;
                long edge_right = nof_h_edges + i * nodes_per_row + (j + 1);
                long edge_diag = nof_h_edges + nof_v_edges + i * n + j;

                // For elements: store the process number in the affiliation field
                // For boundary and fixed_nodes edges: leave the affiliation empty
                long process = i * n + j;
                long empty_aff = 0;

                // Write elements
                long elem_index = (i * n + j) * 2;
                elements(elem_index) = Element{node_bl1, node_br2, node_tr3,
                                               edge_bottom, edge_right, edge_diag, process};
                elements(elem_index + 1) = Element{node_bl1, node_tr3, node_tl4,
                                                   edge_diag, edge_top, edge_left, process};

                // Write boundary and boundary edges
                if (i == 0) {
                    long bdry_index = j;
                    boundary(bdry_index) = BoundaryEdge{node_bl1, node_br2,
                                                        edge_bottom, empty_aff};
                    edges(edge_bottom) = Edge{node_bl1, node_br2};
                }
                if (i == m - 1) {
                    long bdry_index = n + 2 * m + j;
                    boundary(bdry_index) = BoundaryEdge{node_tr3, node_tl4,
                                                        edge_top, empty_aff};
                }
                if (j == 0) {
                    long bdry_index = n + i * 2;
                    boundary(bdry_index) = BoundaryEdge{node_bl1, node_tl4,
                                                        edge_left, empty_aff};
                    edges(edge_left) = Edge{node_bl1, node_tl4};
                }
                if (j == n - 1) {
                    long bdry_index = n + i * 2 + 1;
                    boundary(bdry_index) = BoundaryEdge{node_br2, node_tr3,
                                                        edge_right, empty_aff};
                }

                // Write edges
                edges(edge_top) = Edge{node_tr3, node_tl4};
                edges(edge_right) = Edge{node_br2, node_tr3};
                edges(edge_diag) = Edge{node_bl1, node_tr3};
            }
        }
    }
}
