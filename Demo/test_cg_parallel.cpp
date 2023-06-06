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
    int refine_factor = 4;

    Mesh::GlobalMesh global_mesh;
    Mesh::LocalMesh local_mesh;

    if (rank == 0) {
        global_mesh = Mesh::GlobalMesh(m, n);
        global_mesh.Create();
        global_mesh.Refine(refine_factor);
    }
    Skeleton::Skeleton skeleton(m, n, refine_factor, MPI_COMM_WORLD, rank);
    Util::BlasVector global_sol(global_mesh.get_n_nodes());
    if (rank == 0) {
        skeleton.Create(global_mesh);
        Util::SedMatrix global_stiff = global_mesh.CreateStiffness();
        Util::BlasVector global_rhs = global_mesh.CreateRhs(F_vol, g_Neu);
        global_mesh.AddDirichlet(global_stiff, global_rhs, u_D);
        global_sol = Solver::SolveJacobi(global_stiff, global_rhs);
    }

    global_mesh.Scatter(local_mesh, skeleton);

    // Solve parallel
    Util::SedMatrix stiffness = local_mesh.CreateStiffness();
    Util::BlasVector rhs = local_mesh.CreateRhs(F_vol, g_Neu);
    local_mesh.AddDirichlet(stiffness, rhs, u_D);

    //Util::BlasVector sol = Solver::SolveCGParallel(stiffness, rhs, skeleton);
    Util::BlasVector sol = Solver::SolveJacobiParallel(stiffness, rhs, skeleton);

    Util::BlasVector local_sol_accum;
    skeleton.GatherAccumulatedVector(sol, local_sol_accum);

    // Output
    if (rank == 0) {
        global_sol.Axpy(-1, local_sol_accum);
        printf("\n=========== Serial vs Parallel Solution ===========");
        local_sol_accum.Print();
    }

    MPI_Finalize();
}
