#include "hpc.h"
#include "Mesh.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    index brief = 0;

    RectangularMesh coarse_mesh(m, n);
    mesh_print(&coarse_mesh, brief);
}
