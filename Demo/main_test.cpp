#include "hpc.hpp"
#include <cstdio>

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

int main() {
    int m = 1;
    int n = 2;

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
}
