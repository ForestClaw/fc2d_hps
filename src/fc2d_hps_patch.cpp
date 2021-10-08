#include "fc2d_hps_patch.hpp"

fc2d_hps_patch::fc2d_hps_patch() {}

fc2d_hps_patch::fc2d_hps_patch(fc2d_hps_patchgrid patch_grid, int ID, int level, bool is_leaf) :
	ID(ID), level(level), is_leaf(is_leaf), grid(patch_grid)
		{}

