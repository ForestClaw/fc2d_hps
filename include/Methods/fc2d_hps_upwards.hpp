#ifndef FC2D_HPS_UPWARDS_HPP_
#define FC2D_HPS_UPWARDS_HPP_

#include "fc2d_hps_methods.hpp"

// Globals
extern patch_tree quadtree;
extern int current_ID;
extern std::vector<fc2d_hps_matrix<double>> T_cache;

void visit_set_particular_data_leaves(fc2d_hps_patch& patch);
void visit_upwards(fc2d_hps_patch& patch);
void fc2d_hps_upwards(struct fclaw2d_global* glob);

#endif // FC2D_HPS_UPWARDS_HPP_