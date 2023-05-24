#include "hpc.hpp"
#include <cstdio>

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
    skeleton.CreateLocal(1, local2global, length, 1);
    printf("=== GLOBAL:\n");
    skeleton.Print();
}
