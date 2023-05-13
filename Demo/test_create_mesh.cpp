#include "hpc.h"
#include "mesh.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    index brief = 0;

    RectangularMesh coarse_mesh(m, n);
    coarse_mesh.create();
    mesh_print(coarse_mesh, brief);
}
