#include <cstdio>

#include "hpc.hpp"

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
    return ( 5.0 );
}

int main(int argc, char **argv) {

    int m = 12;
    int n = 8;
    int refine_factor = 5;
    
    Mesh::GlobalMesh mesh(m,n);
 
    mesh.Create();
    mesh.Refine(refine_factor);

    Util::SedMatrix stiffness = mesh.CreateStiffness();
    Util::BlasVector rhs = mesh.CreateRhs(F_vol, g_Neu);
    mesh.AddDirichlet(stiffness, rhs, u_D);

    Util::BlasVector sol = Solver::SolveCg(stiffness, rhs);
    //Util::BlasVector sol = Solver::SolveJacobi(stiffness, rhs);

    // Output
    printf("\n=========== Solution ===========");
    sol.Print();

    printf("\nProblem size: %ld\n", mesh.get_n_nodes());
}
