#include "fc2d_hps_split.hpp"

patch_pair split_vertical(fc2d_hps_patch& tau) {

    // Create child patches
    //    Child grids
    fc2d_hps_patchgrid grid_alpha(tau.grid.Nx, tau.grid.Ny/2, tau.grid.x_lower, tau.grid.x_upper, tau.grid.y_lower, (tau.grid.y_lower + tau.grid.y_upper)/2);
    fc2d_hps_patchgrid grid_beta(tau.grid.Nx, tau.grid.Ny/2, tau.grid.x_lower, tau.grid.x_upper, (tau.grid.y_lower + tau.grid.y_upper)/2, tau.grid.y_upper);

    //    Child patches
    fc2d_hps_patch alpha;
    fc2d_hps_patch beta;

    return {alpha, beta};
}

void split_1to4(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

    // Vertical split


    // Horizontal split

}