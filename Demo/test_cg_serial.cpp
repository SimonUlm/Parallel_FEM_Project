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
    
    
    Mesh::GlobalMesh mesh(3,2);
 
    mesh.Create();
    mesh.Refine();

    Util::SedMatrix stiffness = mesh.CreateStiffness();
    Util::BlasVector rhs = mesh.CreateRhs(F_vol, g_Neu);
    mesh.AddDirichlet(stiffness, rhs, u_D);

    Util::BlasVector sol = Solver::SolveCG(stiffness, rhs);
    //Util::BlasVector sol = Solver::SolveJacobi(stiffness, rhs);

    // Output
    printf("\n=========== Right Hand Side ===========");
    rhs.Print();
    printf("\n=========== Solution ===========");
    sol.Print();
}
