#include "hpc.hpp"
#include <cstdio>
#include <iostream>

double F_vol( Mesh::Node& node, long typ )
{
  return ( 0.0 );
}

double g_Neu( Mesh::Node& node, long typ )
{
  return ( node.x * node.y );
}

double u_D( Mesh::Node& node, long typ )
{
  return ( 10.0 );
  // return ( x[0] * x[1] );
}

int main(int argc, char* argv[]) {

    // Problem size
    int m = 1;
    int n = 2;
    
    //int argv_1 = std::stoi(argv[1]);
    //printf("%d\n",argv_1);

    double t1 = Util::get_wall_time();

    // Assemble problem (build mesh and right hand side)
    Mesh::GlobalMesh mesh(m, n);
    mesh.Create();
    Util::SedMatrix stiffness = mesh.CreateStiffness();
    Util::BlasVector rhs = mesh.CreateRhs(F_vol, g_Neu);
    mesh.AddDirichlet(stiffness, rhs, u_D);

    double t2 = Util::get_wall_time();
    double t_assemble = t2 - t1;

    // Print mesh data, stiffness matrix and rhs
    mesh.Print();
    printf("\n=========== RHS ===========");
    rhs.Print();
    Util::GeMatrix test_ge(stiffness);
    test_ge.Print();
    
    // Solve problem
    t1 = Util::get_wall_time();
    Util::BlasVector sol = Solver::SolveCg(stiffness, rhs);
    t2 = Util::get_wall_time();
    double t_sol = t2 - t1;

    printf("\n=========== SOL ===========");
    sol.Print();

    // Compute A*u, compare to rhs
    Util::BlasVector ref_rhs(mesh.n_nodes());

    // ref_rhs = A * u with general matrix vector product
    stiffness.SymSpmv(1.0, sol, 0.0, ref_rhs);

    // Compute error max(abs(b - A * u))
    ref_rhs.Axpy(-1, rhs);
    double max_err = ref_rhs.Amax();

    printf("\n========== ERROR ==========\n");
    printf("err = %.3e\n", max_err);
    printf("\n====== ELAPSED TIME =======\n");
    printf("duration assembling : %f ms\n", t_assemble);
    printf("duration solving    : %f ms\n", t_sol);
}
