#include "hpc.hpp"
#include <cstdio>

double F_vol( double x[2], index typ )
{
  return ( 0.0 );
}

double g_Neu( double x[2], index typ )
{
  return ( x[0] * x[1] );
}

double u_D( double x[2])
{
  return ( 1.0 );
  // return ( x[0] * x[1] );
}

int main() {
    int m = 1;
    int n = 1;

    Mesh::GlobalMesh mesh(m, n);
    mesh.Create();
    mesh.Print();
    
    Util::SedMatrix stiffness = mesh.CreateStiffness();
    
    Util::BlasVector rhs = mesh.CreateRhs(fvol, fNeu);
    
    mesh.AddDirichlet(stiffness, rhs, u_D);
    
    b.Print();
    
}
