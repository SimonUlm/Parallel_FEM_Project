#include <cstdio>

#include "Include/hpc.hpp"
#include "Include/mesh.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    index brief = 0;

    RectangularMesh mesh(m, n);
    mesh.create();
    mesh.refine();
    mesh_print(mesh, brief);
}
