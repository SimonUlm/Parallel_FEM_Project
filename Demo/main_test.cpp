#include "hpc.hpp"

int main(int argc, char **argv) {
    
    
    Mesh::GlobalMesh mesh(2,3);
 
    mesh.Create();
    mesh.Refine();
    
    
    Util::SedMatrix stiffness_sed = mesh.CreateStiffness();
    Util::BlasVector b(mesh.get_n_nodes());
    b.Init();
    Util::BlasVector sol(mesh.get_n_nodes());

    stiffness_sed.SolveCg(b, sol);

    // Output
    b.Print();
    Util::GeMatrix stiffness_ge(stiffness_sed, true);
    stiffness_ge.Print();
    sol.Print();
}
