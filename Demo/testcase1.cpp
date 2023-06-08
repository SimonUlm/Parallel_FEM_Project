#include <cstdio>
#include <iostream>
#include <unistd.h>

#include "hpc.hpp"

#ifndef MIN_REFINES
#define MIN_REFINES 0
#endif

double F_vol(Mesh::Node &node, long typ) {
    return (3.0);
}

double g_Neu(Mesh::Node &node, long typ) {
    return (node.x * node.y);
}

double u_D(Mesh::Node &node, long typ) {
    return (1.0);
    // return ( x[0] * x[1] );
}

int main(int argc, char *argv[]) {

    // Problem size
    int m = 6;
    int n = 8;
    int refine_factor = std::stoi(argv[1]);

    if (refine_factor == MIN_REFINES) {
        printf("Problem size : %dx%d grid\n", m, n);
        printf("Test case    : CG serial\n");
        printf("------------------------------------------\n\n");
    }

    double error_out;

    // Time stamps
    double t_setup1, t_setup2, t_setup;
    double t_solve1, t_solve2, t_solve;
    double t_assemble1, t_assemble2, t_assemble;
    double t_total;

    Mesh::GlobalMesh global_mesh;

    printf("----------\n");
    printf("refine = %d\n", refine_factor);
    t_setup1 = Util::get_wall_time();
    global_mesh = Mesh::GlobalMesh(m, n);
    global_mesh.Create();
    global_mesh.Refine(refine_factor);
    t_setup2 = Util::get_wall_time();
    t_setup = t_setup2 - t_setup1;
    printf("DOF    = %ld\n", global_mesh.n_nodes());
    printf("----------\n\n");
    printf("t_setup_problem = %f s\n", t_setup);

    t_assemble1 = Util::get_wall_time();
    Util::SedMatrix stiffness = global_mesh.CreateStiffness();
    Util::BlasVector rhs = global_mesh.CreateRhs(F_vol, g_Neu);
    global_mesh.AddDirichlet(stiffness, rhs, u_D);
    t_assemble2 = Util::get_wall_time();
    t_assemble = t_assemble2 - t_assemble1;
    printf("t_assemble      = %f s\n", t_assemble);

    t_solve1 = Util::get_wall_time();
    Util::BlasVector sol = Solver::SolveCg(stiffness, rhs, error_out);
    t_solve2 = Util::get_wall_time();
    t_solve = t_solve2 - t_solve1;
    printf("t_solve         = %f s\n", t_solve);

    t_total = t_setup + t_assemble + t_solve;
    printf("t_total         = %f s\n", t_total);
    printf("----------\n");


    printf("Global residual (norm) = %.3e\n", error_out);
    printf("=======================================\n\n\n");
}
