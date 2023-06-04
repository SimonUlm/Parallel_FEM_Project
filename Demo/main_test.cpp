#include "hpc.hpp"

int main(int argc, char **argv) {
    
    
    Mesh::GlobalMesh mesh(3,4);
 
    mesh.Create();
    mesh.Refine();
    
    
    Util::SedMatrix sed = mesh.CreateStiffness();
    
    sed.Print();
    
    Util::GeMatrix gematrix(sed, true);
    
    gematrix.Print();
    

}
