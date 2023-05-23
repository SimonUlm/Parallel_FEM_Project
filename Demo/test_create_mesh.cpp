#include "hpc.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    GlobalMesh mesh(m, n);
    mesh.Create(Node{1, 1}, Node{2 ,2});
    mesh.Refine();
    mesh.Refine();
    Skeleton skeleton(m, n, mesh.get_refine_factor());
    skeleton.Create(mesh);
    long* local2global = 0;
    long length = 0;
    skeleton.CreateLocal(1, local2global, length);
    skeleton.Print();
}
