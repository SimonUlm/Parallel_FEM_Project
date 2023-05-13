#ifndef HPC2_MESH_HPP
#define HPC2_MESH_HPP

#include <cassert>

namespace Mesh{

    class RectangularMesh {
        // Declare data types
    private:
        struct Node {
            double x;
            double y;
        };

        struct Element {
            long n1;
            long n2;
            long n3;
            long m1;
            long m2;
            long m3;
            long t;
        };

        struct Edge {
            long n1;
            long n2;
        };

        struct BoundaryEdge {
            long n1;
            long n2;
            long m;
            long t;
        };

        template<typename T>
        class List {
        public:
            long count;

            explicit List(long count) :
                    count(count), data(new T[count]) {
            }
            ~List() {
                delete [] data;
            }
            List(const List &) = delete;

            T & operator()(long index) {
                assert(index < count);
                return data[index];
            }
        private:
            T *data;
        };

        // Defines actual class
    public:
        List<Node> nodes;
        List<Element> elements;
        List<Edge> edges;
        List<BoundaryEdge> boundary;
        List<long> fixed_nodes;

        // Support old mesh struct:
        ptrdiff_t ncoord;
        ptrdiff_t nelem;
        ptrdiff_t nedges;
        ptrdiff_t nbdry;
        ptrdiff_t nfixed;
        double *coord;
        ptrdiff_t *elem;
        ptrdiff_t *bdry;
        ptrdiff_t *edge2no;
        ptrdiff_t *fixed;

        RectangularMesh(long m, long n) :
                nodes((m + 1) * (n + 1)), elements(2 * m * n), edges(3 * m * n + m + n),
                boundary(2 * (m + n)), fixed_nodes(m + n + 1) {
            create(m, n);

            // Convert into old struct:
            ncoord = (ptrdiff_t) nodes.count;
            nelem = (ptrdiff_t) elements.count;
            nedges = (ptrdiff_t) edges.count;
            nbdry = (ptrdiff_t) boundary.count;
            nfixed = (ptrdiff_t) fixed_nodes.count;
            coord =  &nodes(0).x;
            elem = (ptrdiff_t *) &elements(0).n1;
            bdry = (ptrdiff_t *) &boundary(0).n1;
            edge2no = (ptrdiff_t *) &edges(0).n1;
            fixed = (ptrdiff_t *) &fixed_nodes(0);
        }

    private:
        void create(long m, long n)
        {
            long node_index = 0;
            long fixed_index = 0;

            long nof_h_edges = (m + 1) * n; // number of horizontal edges
            long nof_v_edges = m * (n + 1); // number of vertical edges
            long nodes_per_row = n + 1;

            // Iterating over rectangles (2 elements at a time)
            for (long i = 0; i < m; ++i) {
                for (long j = 0; j < n; ++j) {

                    // Coordinates of bottom left node
                    double x = 0;
                    double y = 0;
                    if (i == 0 && j == 0) {
                        // Write bottom left node of rectangle for root rectangle
                        nodes(node_index).x = x;
                        nodes(node_index).y = y;
                        fixed_nodes(fixed_index++) = node_index;
                    }

                    if (j == 0) {
                        // Write top left node for first column
                        node_index = nodes_per_row * (i+1);
                        nodes(node_index).x = x;
                        nodes(node_index).y = (double) (i+1) / m;
                        fixed_nodes(fixed_index++) = node_index;
                    }

                    if (i == 0) {
                        // Write bottom right node for first row
                        node_index = j + 1;
                        nodes(node_index).x = (double) (j+1) / n;
                        nodes(node_index).y = y;
                        fixed_nodes(fixed_index++) = node_index;
                    }

                    // Write upper right node for each rectangle
                    x = (double) (j+1) / n;
                    y = (double) (i+1) / m;
                    node_index = nodes_per_row * (i+1) + (j+1);
                    nodes(node_index).x = x;
                    nodes(node_index).y = y;


                    // Four nodes numbers
                    long node_bl1 = i * nodes_per_row + j;            // bottom left node
                    long node_br2 = i * nodes_per_row + j + 1;        // bottom right node
                    long node_tr3 = (i+1) * nodes_per_row + j + 1;    // top right node
                    long node_tl4 = (i+1) * nodes_per_row + j;        // top left node

                    // Five edge numbers
                    long edge_bottom = i * n + j;
                    long edge_top = (i+1) * n + j;
                    long edge_left = nof_h_edges + i * nodes_per_row + j;
                    long edge_right = nof_h_edges + i * nodes_per_row + (j+1);
                    long edge_diag = nof_h_edges + nof_v_edges + i * n + j;

                    // For elements: store the process number in the affiliation field
                    // For boundary and fixed_nodes edges: leave the affiliation empty
                    long process = i * n + j;
                    long empty_aff = 0;

                    // Write elements
                    long elem_index = (i * n + j) * 2;
                    write_elem(elem_index, node_bl1, node_br2, node_tr3,
                               edge_bottom, edge_right, edge_diag, process);
                    write_elem(elem_index + 1, node_bl1, node_tr3, node_tl4,
                               edge_diag, edge_top, edge_left, process);

                    // Write boundary and boundary edges
                    if (i == 0) {
                        long bdry_index = j;
                        write_bdry(bdry_index, node_bl1, node_br2, edge_bottom, empty_aff);
                        write_edge(edge_bottom, node_bl1, node_br2);
                    }
                    if (i == m-1) {
                        long bdry_index = n + 2 * m + j;
                        write_bdry(bdry_index, node_tr3, node_tl4, edge_top, empty_aff);
                    }
                    if (j == 0) {
                        long bdry_index = n + i * 2;
                        write_bdry(bdry_index, node_bl1, node_tl4, edge_left, empty_aff);
                        write_edge(edge_left, node_bl1, node_tl4);
                    }
                    if (j == n-1) {
                        long bdry_index = n + i * 2 + 1;
                        write_bdry(bdry_index, node_br2, node_tr3, edge_right, empty_aff);
                    }

                    // Write edges
                    write_edge(edge_top, node_tr3, node_tl4);
                    write_edge(edge_right, node_br2, node_tr3);
                    write_edge(edge_diag, node_bl1, node_tr3);
                }
            }
        }

        void write_elem(long elem_index,
                        long n1, long n2, long n3,
                        long m1, long m2, long m3,
                        long t)
        {
            elements(elem_index).n1 = n1;
            elements(elem_index).n2 = n2;
            elements(elem_index).n3 = n3;
            elements(elem_index).m1 = m1;
            elements(elem_index).m2 = m2;
            elements(elem_index).m3 = m3;
            elements(elem_index).t = t;
        }

        void write_bdry(long bdry_index,
                        long n1, long n2,
                        long m, long t)
        {
            boundary(bdry_index).n1 = n1;
            boundary(bdry_index).n2 = n2;
            boundary(bdry_index).m = m;
            boundary(bdry_index).t = t;
        }

        void write_edge(long edge_index, long n1, long n2)
        {
            edges(edge_index).n1 = n1;
            edges(edge_index).n2 = n2;
        }
    };
}

#endif //HPC2_MESH_HPP
