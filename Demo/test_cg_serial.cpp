#include "hpc.hpp"

int main(int argc, char **argv) {
    
    
    Mesh::GlobalMesh mesh(2,3);
 
    mesh.Create();
    //mesh.Refine();

    Util::SedMatrix test_sed(mesh.get_n_nodes(), mesh.get_n_nodes() * mesh.get_n_nodes() - 2);
    test_sed.InitDenseSpd();
    Util::BlasVector b(mesh.get_n_nodes());
    b.Init();

    //test_sed.SolveCg(b, sol);
    Util::BlasVector sol = Solver::SolveCG(test_sed, b);

    // Output
    b.Print();
    Util::GeMatrix test_ge(test_sed);
    test_ge.Print();
    sol.Print();
}
