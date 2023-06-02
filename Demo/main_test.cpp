#include "hpc.hpp"

int main(int argc, char **argv) {
    Mesh::GlobalMesh mesh(1,1);
    
    mesh.Create();
    mesh.Print();
    
    
    Util::SedMatrix sed = mesh.CreateStiffness();
    
    sed.Print();
    
    Util::GeMatrix gematrix(sed, true);
    
    gematrix.Print();
    
    
    
    
}
