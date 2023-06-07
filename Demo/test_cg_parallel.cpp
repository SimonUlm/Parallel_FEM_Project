#include <cstdio>

#include "hpc.hpp"

double F_vol( Mesh::Node& node, long typ )
{
    return ( 0.0 );
}

double g_Neu( Mesh::Node& node, long typ )
{
    return ( node.x * node.y );
}

double u_D( Mesh::Node& node, long typ )
{
    return ( 5.0 );
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    int m = 3;
    int n = 2;
    int refine_factor = 5;

    Mesh::GlobalMesh global_mesh;
    Mesh::LocalMesh local_mesh;
    Skeleton::Skeleton skeleton(m, n, refine_factor, MPI_COMM_WORLD, rank);

    // Create problem
    if (rank == 0) {
        global_mesh = Mesh::GlobalMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine(refine_factor);
        skeleton.Create(global_mesh);
    }

    // Scatter problem
    global_mesh.Scatter(local_mesh, skeleton);

    // Solve parallel
    Util::SedMatrix stiffness = local_mesh.CreateStiffness();
    Util::BlasVector rhs = local_mesh.CreateRhs(F_vol, g_Neu);
    local_mesh.AddDirichlet(stiffness, rhs, u_D);

    //Util::BlasVector sol = Solver::SolveCgParallel(stiffness, rhs, skeleton);
    Util::BlasVector sol = Solver::SolveJacobiParallel(stiffness, rhs, skeleton);

    Util::BlasVector local_sol_accum;
    skeleton.GatherAccumulatedVector(sol, local_sol_accum);

    // Output
    if (rank == 0) {
        printf("\n=========== Parallel Solution ===========");
        local_sol_accum.Print();
    }

    MPI_Finalize();
}
