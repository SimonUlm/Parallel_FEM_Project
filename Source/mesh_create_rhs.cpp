#include <cmath>

#include "hpc.hpp"
#include "mesh_objects.hpp"

namespace Mesh {

    void mesh_Neumann(Node &n1, Node &n2, long typ, double (*fc)(Node &, long), double m[2]) {
        Node center_of_edge = n1 + n2;
        center_of_edge = center_of_edge * 0.5;

        Node n21 = n2 - n1;
        double edge_len = sqrt(n21.x * n21.x + n21.y * n21.y);

        double fac = 0.5 * fc(center_of_edge, typ) * edge_len;

        m[0] = fac;
        m[1] = fac;
    }

    void mesh_vol_elem(Node &n1, Node &n2, Node &n3,
                       long typ, double (*fc)(Node &, long),
                       double m[3]) {
        Node n32 = n3 - n2;
        Node n13 = n1 - n3;
        Node n21 = n2 - n1;

        Node center_of_mass = n1 + n2 + n3;
        center_of_mass = center_of_mass * (1.0 / 3);

        // calculate determinant to approximate integral of function f over triangle elem
        double fac = fc(center_of_mass, typ) * (n13.x * n21.y - n21.x * n13.y) / 6.0;

        // Add values for Neumann condition
        for (int i = 0; i < 3; i++) {
            m[i] = fac;
        }
    }

    Util::BlasVector Mesh::CreateRhs(double (*fvol)(Node &, long),
                                     double (*fNeu)(Node &, long)) {
        // For convenience
        long n_nodes = nodes_.count();
        long n_elem = elements_.count();
        long n_bdry = boundary_.count();

        // Calculate the values for the integral (f,\varphi)_\Omega
        long ind[3];
        double bx[3];
        Util::BlasVector b(n_nodes);
        for (long k = 0; k < n_elem; ++k) {
            // Get element nodes_
            ind[0] = elements_(k).n1;
            ind[1] = elements_(k).n2;
            ind[2] = elements_(k).n3;
            // Calculate the integral over k-th element
            mesh_vol_elem(nodes_(ind[0]), nodes_(ind[1]), nodes_(ind[2]), elements_(k).t, fvol, bx);
            for (long j = 0; j < 3; ++j) {
                b(ind[j]) += bx[j];
            }
        }

        // Boundary integrals (Neumann boundaries)
        for (long k = 0; k < n_bdry; ++k) {
            if (boundary_(k).t == 0) {
                // Get edge nodes_
                ind[0] = boundary_(k).n1;
                ind[1] = boundary_(k).n2;

                // Calculate integral over selected Neumann edge
                mesh_Neumann(nodes_(ind[0]), nodes_(ind[1]), boundary_(k).t, fNeu, bx);

                // Add up the values to global rhs vector
                for (long j = 0; j < 2; ++j) {
                    b(ind[j]) += bx[j];
                }
            }
        }
        return b;
    }

    // returns the vector uD (see serial fem paper), update of right hand side is to be performed manually or in other function
    void Mesh::AddDirichlet(Util::SedMatrix &stiff_matrix,
                            Util::BlasVector &b,
                            double (*fDir)(Node &, long)) {
        // For convenience
        long n_nodes = nodes_.count();
        long n_elem = elements_.count();
        long n_bdry = boundary_.count();
        long n_fxd = fixed_nodes_.count();

        Util::BlasVector u(n_nodes);

        // Boundary integrals (Dirichlet boundaries)
        long ind;

        for (long k = 0; k < n_fxd; ++k) {
            // Check if k-th boundary_ is part of dirichlet boundaries
            ind = fixed_nodes_(k);
            u(ind) = fDir(nodes_(ind), 1);
        }

        // b <- b - A * u,  first modification of rhs
        stiff_matrix.SymSpmv(-1, u, 1, b);

        // Adjust stiffness matrix
        for (long k = 0; k < n_fxd; ++k) {
            ind = fixed_nodes_(k);

            // set diagonal to 1
            stiff_matrix(ind) = 1;

            // set column to zero
            stiff_matrix.zero_col(ind);

            // set row to zero
            stiff_matrix.zero_row(ind);
        }

        // Second modification of rhs
        for (long k = 0; k < n_fxd; ++k) {
            // Check if k-th boundary_ is part of dirichlet boundaries
            ind = fixed_nodes_(k);
            b(ind) = u(ind);
        }


    }
}
