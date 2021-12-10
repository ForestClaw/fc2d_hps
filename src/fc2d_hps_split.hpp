#ifndef FC2D_HPS_SPLIT_HPP_
#define FC2D_HPS_SPLIT_HPP_
#pragma once

#include <utility>
#include "fc2d_hps_patch.hpp"

typedef std::pair<fc2d_hps_patch, fc2d_hps_patch> patch_pair;

patch_pair split_vertical(fc2d_hps_patch& tau);
patch_pair split_horizonal(fc2d_hps_patch& tau);
void split_1to4(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3);

#endif // FC2D_HPS_SPLIT_HPP_