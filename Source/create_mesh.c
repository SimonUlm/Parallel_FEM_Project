#include "hpc.h"

// Declare help functions
void write_elem(mesh *wMesh, index elem,
                index n1, index n2, index n3,
                index m1, index m2, index m3,
                index t1);
void write_bdry(mesh *wMesh, index bdry,
                index n1, index n2,
                index m1, index t1);
void write_edge(mesh *wMesh, index edge, index n1, index n2);
void write_fixed(mesh *wMesh, index n);



// Main
mesh *create_rect_mesh(index m, index n) 
{
    index ncoord = (m + 1) * (n + 1);
    index nelem = 2 * m * n;
    index nedges = 3 * m * n + m + n;
    index nbdry = 2 * (m + n);
    index nfixed = m + n + 1;

    mesh *new_mesh = mesh_alloc_with_edges(ncoord, nelem, nbdry, nedges, nfixed);

    index node_index = 0;

    index nof_h_edges = (m + 1) * n; // number of horizontal edges
    index nof_v_edges = m * (n + 1); // number of vertical edges
    index nodes_per_row = n + 1;

    // Iterating over rectangles (2 elements at a time)
    for (index i = 0; i < m; ++i) {
        for (index j = 0; j < n; ++j) {

            // Coordinates of bottom left node
            double x = 0;
            double y = 0;
            if (i == 0 && j == 0) {
            	// Write bottom left node of rectangle for root rectangle
                node_index = 0;
            	new_mesh->coord[node_index * 2] = x;
            	new_mesh->coord[node_index * 2 + 1] = y;
                write_fixed(new_mesh, node_index);
            }
            
            if (j == 0) {
            	// Write top left node for first column
                node_index = nodes_per_row * (i+1);
            	new_mesh->coord[node_index * 2] = x;
            	new_mesh->coord[node_index * 2 + 1] = (double) (i+1) / m;
                write_fixed(new_mesh, node_index);
            }
            
            if (i == 0) {
            	// Write bottom right node for first row
                node_index = j + 1;
            	new_mesh->coord[node_index * 2] = (double) (j+1) / n;
            	new_mesh->coord[node_index * 2 + 1] = y;
                write_fixed(new_mesh, node_index);
            }
            
            // Write upper right node for each rectangle
            x = (double) (j+1) / n;
            y = (double) (i+1) / m;
            node_index = nodes_per_row * (i+1) + (j+1);
            new_mesh->coord[node_index * 2] = x;
            new_mesh->coord[node_index * 2 + 1] = y;


            // Four nodes numbers
            index node_bl1 = i * nodes_per_row + j;              // bottom left node:
            index node_br2 = i * nodes_per_row + j + 1;          // bottom right node:
            index node_tr3 = (i+1) * nodes_per_row + j + 1;    // top right node:
            index node_tl4 = (i+1) * nodes_per_row + j;        // top left node:

            // Five edge numbers
            index edge_bottom = i * n + j;
            index edge_top = (i+1) * n + j;
            index edge_left = nof_h_edges + i * nodes_per_row + j;
            index edge_right = nof_h_edges + i * nodes_per_row + (j+1);
            index edge_diag = nof_h_edges + nof_v_edges + i * n + j;

            // For elements: store the process number in the empty_aff field
            // For boundaries and fixed edges: leave the empty_aff empty
            index process = i * n + j;
            index empty_aff = 0;

            // Write elements
            index elem_index = (i * n + j) * 2;
            write_elem(new_mesh, elem_index, node_bl1, node_br2, node_tr3,
                       edge_bottom, edge_right, edge_diag, process);
            write_elem(new_mesh, elem_index + 1, node_bl1, node_tr3, node_tl4,
                       edge_diag, edge_top, edge_left, process);

            // Write boundaries and boundary edges
            if (i == 0) { // Dirichlet
                index bdry_index = j;
                write_bdry(new_mesh, bdry_index, node_bl1, node_br2, edge_bottom, empty_aff);
                write_edge(new_mesh, edge_bottom, node_bl1, node_br2);
            }
            if (i == m-1) {
                index bdry_index = n + 2 * m + j;
                write_bdry(new_mesh, bdry_index, node_tr3, node_tl4, edge_top, empty_aff);
            }
            if (j == 0) { // Dirichlet
                index bdry_index = n + i * 2;
                write_bdry(new_mesh, bdry_index, node_bl1, node_tl4, edge_left, empty_aff);
                write_edge(new_mesh, edge_left, node_bl1, node_tl4);
            }
            if (j == n-1) {
                index bdry_index = n + i * 2 + 1;
                write_bdry(new_mesh, bdry_index, node_br2, node_tr3, edge_right, empty_aff);
            }

            // Write edges
            write_edge(new_mesh, edge_top, node_tr3, node_tl4);
            write_edge(new_mesh, edge_right, node_br2, node_tr3);
            write_edge(new_mesh, edge_diag, node_bl1, node_tr3);
        }
    }

    return new_mesh;
}



// Define help functions
void write_elem(mesh *wMesh, index elem,
                index n1, index n2, index n3,
                index m1, index m2, index m3,
                index t1) 
{
    wMesh->elem[elem * 7] = n1;
    wMesh->elem[elem * 7 + 1] = n2;
    wMesh->elem[elem * 7 + 2] = n3;
    wMesh->elem[elem * 7 + 3] = m1;
    wMesh->elem[elem * 7 + 4] = m2;
    wMesh->elem[elem * 7 + 5] = m3;
    wMesh->elem[elem * 7 + 6] = t1;
}

void write_bdry(mesh *wMesh, index bdry,
                index n1, index n2,
                index m1, index t1) 
{
    wMesh->bdry[bdry * 4] = n1;
    wMesh->bdry[bdry * 4 + 1] = n2;
    wMesh->bdry[bdry * 4 + 2] = m1;
    wMesh->bdry[bdry * 4 + 3] = t1;
}

void write_edge(mesh *wMesh, index edge, index n1, index n2) 
{
    wMesh->edge2no[edge * 2] = n1;
    wMesh->edge2no[edge * 2 + 1] = n2;
}

void write_fixed(mesh *wMesh, index n)
{
    static index fixed = 0;
    wMesh->fixed[fixed] = n;
    fixed++;
}


