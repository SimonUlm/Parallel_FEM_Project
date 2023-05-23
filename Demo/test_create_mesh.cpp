#include "hpc.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    GlobalMesh mesh(m, n);
    mesh.Create(Node{1, 1}, Node{2 ,2});
    mesh.Refine();
    mesh.Refine();
    mesh.Print();
}
