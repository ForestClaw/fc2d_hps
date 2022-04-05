#include "Methods/fc2d_hps_build.hpp"

void visit_build_dtn(fc2d_hps_patch& patch) {

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
            if (T_cache[patch.level].rows == 0 && T_cache[patch.level].cols == 0) {
                // T for patch level is not set; build from patch solver
                T_cache[patch.level] = FISHPACK_solver.build_dtn(patch.grid);
            }
            patch.T = T_cache[patch.level];
        }
        else {
            patch.T = FISHPACK_solver.build_dtn(patch.grid);
        }

        // @TODO: Add to upwards pass stage of HPS method
        // if (hps_opt->nonhomogeneous_rhs) {
        // Set particular solution
        //    Get RHS data and set to f
        // fclaw2d_domain_t* domain = glob->domain;
        // fclaw2d_patch_t* fc_patch = &(domain->blocks->patches[patch.ID]);
        // int mfields, meqn;
        // double* rhs;
        // double* q;
        // fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
        // fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
        // patch.f = fc2d_hps_vector<double>(patch.grid.Nx * patch.grid.Ny);
        // for (int i = 0; i < patch.grid.Nx; i++) {
        //     for (int j = 0; j < patch.grid.Ny; j++) {
        //         int idx = j + i*patch.grid.Ny;
        //         int idx_T = i + j*patch.grid.Nx;
        //         patch.f[idx] = rhs[idx_T];
        //     }
        // }

        // //    Compute and set particular solution
        // fc2d_hps_vector<double> g_zero(2*patch.grid.Nx + 2*patch.grid.Ny, 0);

        // // Set Neumann data for particular solution
        // patch.h = FISHPACK_solver.dtn(patch.grid, g_zero, patch.f);
        // // }

    }
    return;
}

void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    // printf("4to1 merge: alpha = %i, beta = %i, gamma = %i, omega = %i\n", alpha.ID, beta.ID, gamma.ID, omega.ID);
    merge_4to1(tau, alpha, beta, gamma, omega);

    fc2d_hps_FISHPACK_solver solver;
    fc2d_hps_matrix<double> T_parent = solver.build_dtn(tau.grid);

    double max_diff = 0;
    for (int i = 0; i < tau.T.rows; i++) {
        for (int j = 0; j < tau.T.cols; j++) {
            double diff = fabs(tau.T(i,j) - T_parent(i,j));
            // printf("%i    %i    %16.8e    %16.8e    %16.8e    %16.8e\n", i, j, tau.T(i,j), T_root(i,j), diff, max_diff);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
    }

    printf("max difference for patch at [%.2f, %.2f] to [%.2f, %.2f] = %16.8e\n", tau.grid.x_lower, tau.grid.y_lower, tau.grid.x_upper, tau.grid.y_upper, max_diff);
}

void fc2d_hps_build(fclaw2d_global_t *glob)
{
    fclaw_global_essentialf("Begin HPS build\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance();

    // Build DtN matrix on leaf patches
    quadtree->traverse(visit_build_dtn);

    // Merge DtN's up tree and build solution
    quadtree->merge(visit_merge);

    fclaw_global_essentialf("End HPS build\n");
}