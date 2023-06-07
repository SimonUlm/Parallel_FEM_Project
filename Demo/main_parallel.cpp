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

    // Initialise MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Define problem size
    int m = 3;
    int n = 2;
    int refine_factor = 5;

    // Define variables
    Mesh::GlobalMesh global_mesh;
    Mesh::LocalMesh local_mesh;
    Skeleton::Skeleton skeleton(m, n, refine_factor, MPI_COMM_WORLD, rank);

    // Create problem (root only)
    if (rank == 0) {
        global_mesh.Create(m, n);
        global_mesh.Refine(refine_factor);
        skeleton.Create(global_mesh);
    }

    // Scatter problem
    global_mesh.Scatter(local_mesh, skeleton);

    // Assemble matrix and rhs
    Util::SedMatrix stiffness = local_mesh.CreateStiffness();
    Util::BlasVector rhs = local_mesh.CreateRhs(F_vol, g_Neu);
    local_mesh.AddDirichlet(stiffness, rhs, u_D);

    // Solve problem
    double error;
    Util::BlasVector local_sol = Solver::SolveCgParallel(stiffness, rhs, skeleton, error);

    // Gather solution
    Util::BlasVector global_sol;
    skeleton.GatherAccumulatedVector(local_sol, global_sol);

    // Output (root only)
    if (rank == 0) {
        printf("\n=========== Solution Vector ===========");
        global_sol.Print();
        printf("Error = %lf\n", error);
    }

    MPI_Finalize();
}