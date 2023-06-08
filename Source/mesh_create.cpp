#include "hpc.hpp"

namespace Mesh {

    /*
     * Creation of an unrefined global mesh as a rectangle which is defined by the
     * two given nodes_ (default (0,0) (1,1)
     *
     * bottom_left_node: node to determine coordinates of bottom left point of mesh
     * top_right_node: node to determine coordinates of top right point of mesh
     *
     */
    void GlobalMesh::Create(Node bottom_left_node, Node top_right_node) {
        long node_index = 0;
        long fixed_index = 0;

        long nof_h_edges = (m_ + 1) * n_; // number of horizontal edges
        long nof_v_edges = m_ * (n_ + 1); // number of vertical edges
        long nodes_per_row = n_ + 1;

        double delta_x = top_right_node.x - bottom_left_node.x;
        double delta_y = top_right_node.y - bottom_left_node.y;

        // Iterating over rectangles (2 elements_ at a time)
        for (long i = 0; i < m_; ++i) {
            for (long j = 0; j < n_; ++j) {

                // Coordinates of bottom left node
                if (i == 0 && j == 0) {
                    // Write bottom left node of rectangle for root rectangle
                    nodes_(node_index) = bottom_left_node;
                    fixed_nodes_(fixed_index++) = node_index;
                }

                if (j == 0) {
                    // Write top left node for first column
                    node_index = nodes_per_row * (i + 1);
                    nodes_(node_index).x = bottom_left_node.x;
                    nodes_(node_index).y = bottom_left_node.y + ((double) (i + 1) / m_) * delta_y;
                    fixed_nodes_(fixed_index++) = node_index;
                }

                if (i == 0) {
                    // Write bottom right node for first row
                    node_index = j + 1;
                    nodes_(node_index).x = bottom_left_node.x + ((double) (j + 1) / n_) * delta_x;
                    nodes_(node_index).y = bottom_left_node.y;
                    fixed_nodes_(fixed_index++) = node_index;
                }

                // Write upper right node for each rectangle
                node_index = nodes_per_row * (i + 1) + (j + 1);
                nodes_(node_index).x = bottom_left_node.x + ((double) (j + 1) / n_) * delta_x;
                nodes_(node_index).y = bottom_left_node.y + ((double) (i + 1) / m_) * delta_y;


                // Four nodes_ numbers
                long node_bl1 = i * nodes_per_row + j;            // bottom left node
                long node_br2 = i * nodes_per_row + j + 1;        // bottom right node
                long node_tr3 = (i + 1) * nodes_per_row + j + 1;    // top right node
                long node_tl4 = (i + 1) * nodes_per_row + j;        // top left node

                // Five edge numbers
                long edge_bottom = i * n_ + j;
                long edge_top = (i + 1) * n_ + j;
                long edge_left = nof_h_edges + i * nodes_per_row + j;
                long edge_right = nof_h_edges + i * nodes_per_row + (j + 1);
                long edge_diag = nof_h_edges + nof_v_edges + i * n_ + j;

                // For elements: store the process number in the affiliation field
                // For boundary: store 1 into the dirichlet boundary edges
                // For boundary and fixed_nodes_ edges_: leave the affiliation empty
                long process = i * n_ + j;
                long dirichlet = 1;
                long empty_aff = 0;


                // Write elements
                long elem_index = (i * n_ + j) * 2;
                elements_(elem_index) = Element{node_bl1, node_br2, node_tr3,
                                                edge_bottom, edge_right, edge_diag, process};
                elements_(elem_index + 1) = Element{node_bl1, node_tr3, node_tl4,
                                                    edge_diag, edge_top, edge_left, process};

                // Write boundary and boundary edges
                if (i == 0) {
                    // boundaries of first processor row are defined as dirichlet
                    long bdry_index = j;
                    boundary_(bdry_index) = BoundaryEdge{node_bl1, node_br2,
                                                         edge_bottom, dirichlet};
                    edges_(edge_bottom) = Edge{node_bl1, node_br2};
                }
                if (i == m_ - 1) {
                    long bdry_index = n_ + 2 * m_ + j;
                    boundary_(bdry_index) = BoundaryEdge{node_tr3, node_tl4,
                                                         edge_top, empty_aff};
                }
                if (j == 0) {
                    // boundaries of first processor col are defined as dirichlet
                    long bdry_index = n_ + i * 2;
                    boundary_(bdry_index) = BoundaryEdge{node_bl1, node_tl4,
                                                         edge_left, dirichlet};
                    edges_(edge_left) = Edge{node_bl1, node_tl4};
                }
                if (j == n_ - 1) {
                    long bdry_index = n_ + i * 2 + 1;
                    boundary_(bdry_index) = BoundaryEdge{node_br2, node_tr3,
                                                         edge_right, empty_aff};
                }

                // Write edges
                edges_(edge_top) = Edge{node_tr3, node_tl4};
                edges_(edge_right) = Edge{node_br2, node_tr3};
                edges_(edge_diag) = Edge{node_bl1, node_tr3};
            }
        }
    }
}
