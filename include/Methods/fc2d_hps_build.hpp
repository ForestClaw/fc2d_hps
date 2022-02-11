#ifndef FC2D_HPS_BUILD_HPP_
#define FC2D_HPS_BUILD_HPP_

#include "fc2d_hps_methods.hpp"

// Global declarations
extern patch_tree quadtree; // use static...?
extern int current_ID;
extern std::vector<fc2d_hps_matrix<double>> T_cache;

static void cb_merge(fclaw2d_global_t *glob, fclaw2d_patch_t *fine_patches, int blockno, int fine0_patchno, void *user);
void visit_leaves(fc2d_hps_patch& patch);
void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega);
void fc2d_hps_build(struct fclaw2d_global* glob);

#endif // FC2D_HPS_BUILD_HPP_