#ifndef FC2D_HPS_MERGE_HPP
#define FC2D_HPS_MERGE_HPP

#include <numeric>
#include <vector>
#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"
#include "fc2d_hps_patch.hpp"

fc2d_hps_patch merge_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
fc2d_hps_patch merge_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
void merge_4to1(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3);

#endif // FC2D_HPS_MERGE_HPP