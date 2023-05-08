#include "hpc.h" 
#include <stdio.h>


int main()
{
    // printf("hello world\n");
    index m = 4;
    index n = 3;
    
    index brief = 0;

    mesh* newMesh = createMesh(m, n);
    printf("Mesh created!\n");
    
    index a = mesh_print(newMesh, brief);
}
