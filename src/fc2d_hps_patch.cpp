#include "fc2d_hps_patch.hpp"

fc2d_hps_patch::fc2d_hps_patch() {}

fc2d_hps_patch::fc2d_hps_patch(fc2d_hps_patchgrid patch_grid, int ID, int level, bool is_leaf) :
	ID(ID), level(level), is_leaf(is_leaf), grid(patch_grid) {
		if (is_leaf) {
			N_patch_side[WEST] = 1;
			N_patch_side[EAST] = 1;
			N_patch_side[SOUTH] = 1;
			N_patch_side[NORTH] = 1;
		}
	}

