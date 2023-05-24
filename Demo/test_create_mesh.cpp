#include "hpc.hpp"
#include <cstdio>

int main() {
    int m = 4;
    int n = 3;

    Mesh::GlobalMesh mesh(m, n);
    mesh.Create();
    mesh.Refine();
    mesh.Refine();
    mesh.Print();
    Skeleton::Skeleton skeleton(m, n, mesh.get_refine_factor());
    skeleton.Create(mesh);
    skeleton.Print();
}
