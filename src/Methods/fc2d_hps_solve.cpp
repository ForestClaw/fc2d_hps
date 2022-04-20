#include <Methods/fc2d_hps_solve.hpp>

// Global declarations
extern std::vector<fc2d_hps_matrix<double>> T_cache;

void visit_split(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    if (tau.is_leaf == false) {
        // printf("splitting...\n");
        split_1to4(tau, alpha, beta, gamma, omega);
    }
}

void visit_patchsolver(fc2d_hps_patch& patch) {
    if (patch.is_leaf == true) {

        // Get options
        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);

        // Get patch solver
        // TODO: Put patch solver in glob...?
        fc2d_hps_FISHPACK_solver patch_solver;

        if (hps_opt->nonhomogeneous_rhs == 0) {
            patch.f = fc2d_hps_vector<double>(patch.grid.Nx * patch.grid.Ny, 0);
        }
        patch.u = patch_solver.solve(patch.grid, patch.g, patch.f);
    }
}

void set_root_boundary_data(fc2d_hps_patch& root_patch) {

    fclaw2d_global_t* glob = (fclaw2d_global_t*) root_patch.user;
    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);
    fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();

    int size_of_g = 2*root_patch.grid.Nx + 2*root_patch.grid.Ny;
    root_patch.g = fc2d_hps_vector<double>(size_of_g, 0);
    fc2d_hps_vector<double> g_west(root_patch.grid.Ny);
    fc2d_hps_vector<double> g_east(root_patch.grid.Ny);
    fc2d_hps_vector<double> g_south(root_patch.grid.Nx);
    fc2d_hps_vector<double> g_north(root_patch.grid.Nx);
    for (int j = 0; j < root_patch.grid.Ny; j++) {
        double y = root_patch.grid.point(YDIM, j);
        g_west[j] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[0], &glob->curr_time, &root_patch.grid.x_lower, &y);
        g_east[j] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[1], &glob->curr_time, &root_patch.grid.x_upper, &y);
    }
    for (int i = 0; i < root_patch.grid.Nx; i++) {
        double x = root_patch.grid.point(XDIM, i);
        g_south[i] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[2], &glob->curr_time, &x, &root_patch.grid.y_lower);
        g_north[i] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[3], &glob->curr_time, &x, &root_patch.grid.y_upper);
    }
    root_patch.g.intract(0*root_patch.grid.Nx, g_west);
    root_patch.g.intract(1*root_patch.grid.Nx, g_east);
    root_patch.g.intract(2*root_patch.grid.Nx, g_south);
    root_patch.g.intract(3*root_patch.grid.Nx, g_north);
    return;

}

void fc2d_hps_solve(fclaw2d_global_t* glob) {
    
    fclaw_global_essentialf("Begin HPS solve\n");

    // Get options
    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);
    fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance();

    // Build Dirichlet data at top level
    fc2d_hps_patch& root_patch = quadtree->data[0];
    set_root_boundary_data(root_patch);

    // Traverse tree from root and apply solution operator or patch solver
    quadtree->split(visit_split);

    // Iterate over leaf nodes and apply patch solver
    quadtree->traverse_preorder(visit_patchsolver);

    fclaw_global_essentialf("End HPS solve\n");
}

void fc2d_hps_fishpack_solve(fclaw2d_global_t* glob) {

    // Get stuff from glob
    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);
    fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();

    FCLAW_ASSERT(fclaw_opt->minlevel == 0 && fclaw_opt->maxlevel = 0);

    // Build grid
    int Nx = clawpatch_opt->mx;
	int Ny = clawpatch_opt->my;
	double x_lower = fclaw_opt->ax;
	double x_upper = fclaw_opt->bx;
	double y_lower = fclaw_opt->ay;
	double y_upper = fclaw_opt->by;
	fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

    // Get Dirichlet data
    int size_of_g = 2*grid.Nx + 2*grid.Ny;
    fc2d_hps_vector<double> g = fc2d_hps_vector<double>(size_of_g, 0);
    fc2d_hps_vector<double> g_west(grid.Ny);
    fc2d_hps_vector<double> g_east(grid.Ny);
    fc2d_hps_vector<double> g_south(grid.Nx);
    fc2d_hps_vector<double> g_north(grid.Nx);
    for (int j = 0; j < grid.Ny; j++) {
        double y = grid.point(YDIM, j);
        g_west[j] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[0], &glob->curr_time, &grid.x_lower, &y);
        g_east[j] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[1], &glob->curr_time, &grid.x_upper, &y);
    }
    for (int i = 0; i < grid.Nx; i++) {
        double x = grid.point(XDIM, i);
        g_south[i] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[2], &glob->curr_time, &x, &grid.y_lower);
        g_north[i] = hps_vt->fort_eval_bc(&hps_opt->boundary_conditions[3], &glob->curr_time, &x, &grid.y_upper);
    }
    g.intract(0*grid.Nx, g_west);
    g.intract(1*grid.Nx, g_east);
    g.intract(2*grid.Nx, g_south);
    g.intract(3*grid.Nx, g_north);

    // Get RHS data
    fclaw2d_patch_t* fc_patch = &(glob->domain->blocks->patches[0]);
    int mfields;
    double* rhs;
    fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
    fc2d_hps_vector<double> f(grid.Nx * grid.Ny);
    for (int i = 0; i < grid.Nx; i++) {
        for (int j = 0; j < grid.Ny; j++) {
            int idx = j + i*grid.Ny;
            int idx_T = i + j*grid.Nx;
            f[idx] = rhs[idx_T];
        }
    }

    // Create patch solver
    fc2d_hps_FISHPACK_solver patch_solver;

    // Solve with FISHPACK
    fc2d_hps_vector<double> u_FISHPACK = patch_solver.solve(grid, g, f);

    // Move into RHS of ForestClaw data
    int meqn;
    double* q;
    fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
    for (int i = 0; i < grid.Nx; i++) {
        for (int j = 0; j < grid.Ny; j++) {
            int idx = j + i*grid.Ny;
            int idx_T = i + j*grid.Nx;
            q[idx_T] = u_FISHPACK[idx];
            rhs[idx_T] = u_FISHPACK[idx];
        }
    }

}