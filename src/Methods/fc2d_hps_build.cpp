#include "Methods/fc2d_hps_build.hpp"

void visit_leaves(fc2d_hps_patch& patch) {

    // TODO: Add level optimizations
    if (patch.is_leaf) {
        // Create patch solver
        fc2d_hps_FISHPACK_solver FISHPACK_solver;

        // Get hps options
        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fc2d_hps_options* hps_opt = fc2d_hps_get_options(glob);

        // Set DtN on leaf
        // @TODO: Implement caching of T for each level
        if (hps_opt->cache_T) {

            // printf("size of cache: %i\n", T_cache.capacity());

            if (T_cache[patch.level].rows == 0 && T_cache[patch.level].cols == 0) {
                // T for patch level is not set; build from patch solver
                T_cache[patch.level] = FISHPACK_solver.build_dtn(patch.grid);
                // printf("size of cache[%i] = %i, %i\n", patch.level, T_cache[patch.level].rows, T_cache[patch.level].cols);
            }
            // patch.T = FISHPACK_solver.build_dtn(patch.grid);
            patch.T = T_cache[patch.level];

            patch.print_info();
        }
        else {
            patch.T = FISHPACK_solver.build_dtn(patch.grid);
        }

        if (hps_opt->nonhomogeneous_rhs) {
            // Set particular solution
            //    Get RHS data and set to f
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

        // patch.print_info();

    }
    return;
}

void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    merge_4to1(tau, alpha, beta, gamma, omega);
}

void fc2d_hps_build(fclaw2d_global_t *glob)
{
    fclaw_global_essentialf("Begin HPS build\n");

    // Build DtN matrix on leaf patches
    quadtree.traverse_inorder(visit_leaves);

    // Merge DtN's up tree and build solution
    quadtree.merge(visit_merge);

    fclaw_global_essentialf("End HPS build\n");
}