#ifndef FC2D_HPS_SPLIT_HPP_
#define FC2D_HPS_SPLIT_HPP_
#pragma once

#include "fc2d_hps_patch.hpp"

void split_vertical(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
void split_horizonal(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
void split_1to4(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3);

#endif // FC2D_HPS_SPLIT_HPP_