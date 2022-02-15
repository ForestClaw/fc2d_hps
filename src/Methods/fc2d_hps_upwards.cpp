#include <Methods/fc2d_hps_upwards.hpp>

void visit_set_particular_data_leaves(fc2d_hps_patch& patch) {

    if (patch.is_leaf) {

        // Create patch solver
        fc2d_hps_FISHPACK_solver FISHPACK_solver;

        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fc2d_hps_options* hps_opt = fc2d_hps_get_options(glob);

        if (hps_opt->nonhomogeneous_rhs) {
            // Set particular solution
            //    Get RHS data and set to patch's f
            fclaw2d_domain_t* domain = glob->domain;
            fclaw2d_patch_t* fc_patch = &(domain->blocks->patches[patch.ID]);
            int mfields;
            double* rhs;
            fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
            patch.f = fc2d_hps_vector<double>(rhs, rhs + patch.grid.Nx*patch.grid.Ny);

            //    Compute and set particular solution
            fc2d_hps_vector<double> g_zero(2*patch.grid.Nx + 2*patch.grid.Ny, 0);
            patch.w = FISHPACK_solver.solve(patch.grid, g_zero, patch.f);

            // Set Neumann data for particular solution
            patch.h = FISHPACK_solver.dtn(patch.grid, g_zero, patch.f);
        }

    }
    // patch.print_info();
}

void fc2d_hps_upwards(fclaw2d_global_t* glob) {
    fclaw_global_essentialf("Begin HPS upwards pass\n");

    // Traverse inorder to set non-homogeneous RHS and BC data into solution vector
    // quadtree.traverse_inorder(visit_upwards);
    fclaw_global_essentialf("End HPS upwards pass\n");
}