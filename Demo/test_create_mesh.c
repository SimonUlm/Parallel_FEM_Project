#include "hpc.h"


int main() {
    index m = 4;
    index n = 3;

    index brief = 0;

    mesh *newMesh = create_rect_mesh(m, n);

    mesh_print(newMesh, brief);
}
