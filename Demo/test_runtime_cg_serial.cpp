    Util::SedMatrix global_stiffness = global_mesh.CreateStiffness();
    Util::BlasVector global_rhs = global_mesh.CreateRhs(F_vol, g_Neu);
    global_mesh.AddDirichlet(global_stiffness, global_rhs, u_D);
    Util::BlasVector global_sol = Solver::SolveCG(global_stiffness, global_rhs);
