#include "hpc.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    GlobalMesh mesh(m, n);
    mesh.Create(Node{1, 1}, Node{2 ,2});
    mesh.Print();
    mesh.Refine();
    long refine_factor = 1;
    Skeleton skeleton(m, n, refine_factor);
    skeleton.Create(mesh);
    skeleton.Print();
}
