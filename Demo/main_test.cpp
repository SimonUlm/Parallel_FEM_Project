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
    int m = 1;
    int n = 2;
    
    int argv_1 = std::stoi(argv[1]);
    printf("%d\n",argv_1);
    Mesh::GlobalMesh mesh(m, n);

    mesh.Create();

    mesh.Print();
    
    Util::SedMatrix stiffness = mesh.CreateStiffness();
    
    Util::BlasVector rhs = mesh.CreateRhs(F_vol, g_Neu);
    
    mesh.AddDirichlet(stiffness, rhs, u_D);
    printf("\n=========== RHS ===========");
    rhs.Print();

    Util::GeMatrix test_ge(stiffness);
    test_ge.Print();
    Util::BlasVector sol = Solver::SolveCG(stiffness, rhs);
    printf("\n=========== SOL ===========");
    sol.Print();

    // Compute A*u, compare to rhs
    Util::BlasVector ref_rhs(mesh.get_n_nodes());
    // ref_rhs = A * u
    stiffness.SymSpmv(1.0, sol, 0.0, ref_rhs);
    // infinity norm of ref_rhs - rhs

    ref_rhs.Axpy(-1, rhs);
    double max_err = ref_rhs.Amax();

    printf("\n=========== ERR ===========\n");
    printf("err = %.3e\n", max_err);
}
