#include <cstdio>

#include "Include/mesh.hpp"

using namespace Mesh;

int main() {
    int m = 4;
    int n = 3;

    RectangularMesh mesh(m, n);
    mesh.create();
    mesh.refine();
    mesh.print();
}
