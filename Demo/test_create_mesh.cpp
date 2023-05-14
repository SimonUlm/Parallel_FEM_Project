#include <cstdio>

#include "hpc.h"
#include "mesh.hpp"

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
