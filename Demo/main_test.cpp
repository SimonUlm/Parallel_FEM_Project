#include "hpc.hpp"

int main(int argc, char **argv) {
    Mesh::GloabalMesh mesh(4,3);
    
    mesh.Create();
    
    
    Util::SedMatrix sed = mesh.CreateStiffness();
    sed.InitOne();
    
    sed.Print();
    
    
}
