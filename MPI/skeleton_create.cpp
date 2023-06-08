#include <cmath>

#include "hpc.hpp"

#ifdef _MPI

namespace Skeleton {

    long mn_r(long m, long r);

    /*
     * Create Skeleton from given mesh
     *
     * Creates the global Skeleton from given Mesh. With the number of processes of the Mesh
     * the cross points and Processes of the ComBorder's are calculated. With the refinement
     * factor the ComBorderNodes are calculated and saved.
     *
     * mesh: Mesh for which the Skeleton should be created.
     *
     */
    void Skeleton::Create(Mesh::GlobalMesh &mesh) {
        long border_index = 0;
        long m = mesh.m();
        long m_n = m + 1;                            // Number of nodes in direction of m
        long n = mesh.n();
        long n_n = n + 1;                            // Number of nodes in direction of n

        long nof_h_edges = (m + 1) * n;              // Number of horizontal edges
        long nof_v_edges = m * (n + 1);              // Number of vertical edges
        long nodes_per_row = n + 1;
        long refine = mesh.refine_factor();          // Number of performed refinements

        // Iterate through processes
        // Initialize always left vertical border and lower horizontal border as new couple
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j < n; ++j) {
                // Skip first process in row
                if (j != 0) {
                    // initialize left vertical couple
                    long c1 = i * nodes_per_row + j;        // lower left Node
                    long c2 = nodes_per_row * (i + 1) + j;    // upper left Node
                    long l_proc = i * n + (j - 1);            // left process
                    long r_proc = i * n + j;                // this process
                    long color = (j - 1) % 2;                    // color

                    com_borders_(border_index).set_entries(
                            border_index, c1, c2, l_proc, r_proc, color
                    );

                    long index = 0;
                    // Calculated new nodes_ for each performed refinement
                    // Check documentation for additional information about the formula
                    for (long r = refine; r > 0; --r) {
                        long node = mn_r(m_n, r - 1) * mn_r(n_n, r - 1) + (mn_r(n_n, r - 1) - 1) * m_n;
                        node += pow(2, r - 1) * j;            // accounting column
                        node += pow(2, r - 1) * n_n * i;        // accounting row
                        for (long k = 0; k < pow(2, r - 1); ++k) {
                            com_border_nodes_.set_entry(border_index, index, node);
                            node += 1;
                            index += 1;
                        }
                    }
                    border_index++;
                }

                // Skip first row of processes
                if (i != 0) {
                    // initialize lower horizontal couple
                    long c1 = nodes_per_row * i + j + 1;    // lower right node
                    long c2 = nodes_per_row * i + j;        // lower left node
                    long l_proc = (i - 1) * n + j;            // lower process
                    long r_proc = i * n + j;                // this process as upper process
                    long color = (i - 1) % 2 + 2;                // color_

                    com_borders_(border_index).set_entries(
                            border_index, c1, c2, l_proc, r_proc, color
                    );

                    long index = 0;
                    // Calculated new nodes_ for each performed refinement
                    // Check documentation for additional information about the formula
                    for (long r = refine; r > 0; --r) {
                        long node = mn_r(m_n, r - 1) * mn_r(n_n, r - 1);
                        node += pow(2, r - 1) * j;            // accounting column
                        node += pow(2, r - 1) * n * i;        // accounting row
                        for (long k = 0; k < pow(2, r - 1); ++k) {
                            com_border_nodes_.set_entry(border_index, index, node);
                            node += 1;
                            index += 1;
                        }
                    }
                    border_index++;
                }
            }
        }
    }

    /*
     * Calculate number of nodes in one dimension after r refinements
     *
     * m: Number of processes in on dimension
     * r: number of performed refinements
     * return: Returns number nodes_ in direction of given dimension
     *
     */
    long mn_r(long m, long r) {
        // m_ = number of nodes_ in one dimension
        return (m - 1) * pow(2, r) + 1;
    }
}

#endif // _MPI
