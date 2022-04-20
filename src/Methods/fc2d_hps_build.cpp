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
    }
    return;
}

void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    merge_4to1(tau, alpha, beta, gamma, omega);
}

void fc2d_hps_build(fclaw2d_global_t *glob)
{
    fclaw_global_essentialf("Begin HPS build\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance();

    // Build DtN matrix on leaf patches
    quadtree->traverse_postorder(visit_build_dtn);

    // Merge DtN's up tree and build solution
    quadtree->merge(visit_merge);

    fclaw_global_essentialf("End HPS build\n");
}