#include <cstdio>
#include <iostream>
#include <mpi.h>
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
    double t_mpi1, t_mpi2, t_mpi;
    t_mpi1 = Util::get_wall_time();
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    if (rank == 0) {
        t_mpi2 = Util::get_wall_time();
        t_mpi = t_mpi2 - t_mpi1;

    }

    double error;

    // Time stamps
    double t_setup1, t_setup2, t_setup;
    double t_scatter1, t_scatter2, t_scatter;
    double t_assemble1, t_assemble2, t_assemble;
    double t_solve1, t_solve2, t_solve;
    double t_gather_sol1, t_gather_sol2, t_gather_sol;
    double t_total;

    // Problem size
    int m = 2;
    int n = 3;
    int refine_factor = std::stoi(argv[1]);

    if (refine_factor == MIN_REFINES && rank == 0) {
        printf("Problem size : %dx%d grid, %d processes\n", m, n, m * n);
        printf("Test case    : CG serial parallel\n");
        printf("------------------------------------------\n\n");
    }

    Mesh::GlobalMesh global_mesh;
    Mesh::LocalMesh local_mesh;
    Skeleton::Skeleton skeleton;

    if (rank == 0) {
        printf("----------\n");
        printf("refine = %d\n", refine_factor);
        t_setup1 = Util::get_wall_time();
        global_mesh = Mesh::GlobalMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine(refine_factor);
        skeleton.Create(global_mesh, m, n);
        t_setup2 = Util::get_wall_time();
        t_setup = t_setup2 - t_setup1;
        printf("DOF    = %ld\n", global_mesh.n_nodes());
        printf("----------\n");
        printf("time MPI Setup  = %f s\n", t_mpi);
        printf("t_setup_problem = %f s\n", t_setup);
    }

    if (rank == 0) { t_scatter1 = Util::get_wall_time(); }
    global_mesh.Scatter(local_mesh, skeleton, MPI_COMM_WORLD, rank);
    if (rank == 0) {
        t_scatter2 = Util::get_wall_time();
        t_scatter = t_scatter2 - t_scatter1;
        printf("t_scatter       = %f s\n", t_scatter);
    }

    // Solve parallel
    if (rank == 0) { t_assemble1 = Util::get_wall_time(); }
    Util::SedMatrix local_stiffness = local_mesh.CreateStiffness();
    Util::BlasVector local_rhs = local_mesh.CreateRhs(F_vol, g_Neu);
    local_mesh.AddDirichlet(local_stiffness, local_rhs, u_D);
    if (rank == 0) {
        t_assemble2 = Util::get_wall_time();
        t_assemble = t_assemble2 - t_assemble1;
        printf("t_assemble      = %f s\n", t_assemble);
    }
    if (rank == 0) { t_solve1 = Util::get_wall_time(); }
    Util::BlasVector sol = Solver::SolveCgParallel(local_stiffness, local_rhs, skeleton, error);
    if (rank == 0) {
        t_solve2 = Util::get_wall_time();
        t_solve = t_solve2 - t_solve1;
        printf("t_solve         = %f s\n", t_solve);
    }

    if (rank == 0) { t_gather_sol1 = Util::get_wall_time(); }
    Util::BlasVector local_sol_accum;
    skeleton.GatherAccumulatedVector(sol, local_sol_accum);
    if (rank == 0) {
        t_gather_sol2 = Util::get_wall_time();
        t_gather_sol = t_gather_sol2 - t_gather_sol1;
        printf("t_gather_sol    = %f s\n", t_gather_sol2 - t_gather_sol1);
        t_total = t_mpi + t_setup + t_scatter + t_assemble + t_solve + t_gather_sol;
        printf("t_total         = %f s\n", t_total);
        printf("----------\n");
        printf("Global residual (norm) = %.3e\n", error);
        printf("=======================================\n\n\n");
    }

    MPI_Finalize();
}
