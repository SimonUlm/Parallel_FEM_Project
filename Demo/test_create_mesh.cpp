#include "hpc.hpp"

int main() {
    int m = 4;
    int n = 3;

    Mesh::GlobalMesh mesh(m, n);
    mesh.Create();
    mesh.Refine();
    mesh.Refine();
    Skeleton::Skeleton skeleton(m, n, mesh.get_refine_factor());
    skeleton.Create(mesh);
    long* local2global = 0;
    long length = 0;
    skeleton.CreateLocal(1, local2global, length);
    skeleton.Print();
}
