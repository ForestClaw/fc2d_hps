#ifndef FC2D_HPS_SOLVE_HPP_
#define FC2D_HPS_SOLVE_HPP_

#include "fc2d_hps_methods.hpp"

extern patch_tree quadtree;
extern int current_ID;

void visit_split(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega);
void visit_patchsolver(fc2d_hps_patch& patch);
void set_root_boundary_data(fc2d_hps_patch& root_patch);
void fc2d_hps_solve(struct fclaw2d_global* glob);

void fc2d_hps_fishpack_solve(fclaw2d_global_t* glob);

#endif /* !FC2D_HPS_SOLVE_HPP_ */

