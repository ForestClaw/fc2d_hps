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

void fc2d_hps_patch::print_info() {
	printf("---------- PATCH INFO ----------\n");
	printf("  ID           = %i\n", this->ID);
	printf("  level        = %i\n", this->level);
	printf("  is_leaf      = %i\n", this->is_leaf);
	printf("  N_cells_leaf = %i\n", this->N_cells_leaf);
	printf("  N_patch_side:\n");
	printf("    WEST  = %i\n", this->N_patch_side[WEST]);
	printf("    EAST  = %i\n", this->N_patch_side[EAST]);
	printf("    SOUTH = %i\n", this->N_patch_side[SOUTH]);
	printf("    NORTH = %i\n", this->N_patch_side[NORTH]);
	printf("  Grid:\n");
	printf("    grid.Nx      = %i\n", this->grid.Nx);
	printf("    grid.Ny      = %i\n", this->grid.Ny);
	printf("    grid.dx      = %11.4e\n", this->grid.dx);
	printf("    grid.dy      = %11.4e\n", this->grid.dy);
	printf("    grid.x_lower = %11.4e\n", this->grid.x_lower);
	printf("    grid.x_upper = %11.4e\n", this->grid.x_upper);
	printf("    grid.y_lower = %11.4e\n", this->grid.y_lower);
	printf("    grid.y_upper = %11.4e\n", this->grid.y_upper);
	printf("  Data Sizes:\n");
	printf("    T       : [%8i, %8i]\n", this->T.rows, this->T.cols);
	printf("    S       : [%8i, %8i]\n", this->S.rows, this->S.cols);
	printf("    S_prime : [%8i, %8i]\n", this->S_prime.rows, this->S_prime.cols);
	printf("    X       : [%8i, %8i]\n", this->X.rows, this->X.cols);
	printf("    H       : [%8i, %8i]\n", this->H.rows, this->H.cols);
	printf("    u       : [%8i]\n", this->u.size());
	printf("    f       : [%8i]\n", this->f.size());
	printf("    g       : [%8i]\n", this->g.size());
	printf("    h       : [%8i]\n", this->h.size());
	printf("    w       : [%8i]\n", this->w.size());
	printf("--------------------------------\n");
}
