#include "Methods/fc2d_hps_setup.hpp"

// patch_tree quadtree;
int current_ID;
std::vector<fc2d_hps_matrix<double>> T_cache(25); // TODO: Create cache class

fc2d_hps_patch init_fn(fc2d_hps_quadtree<fc2d_hps_patch>* quadtree, int level, int idx, void* user) {

	// Get glob from user and get options
	fclaw2d_global_t* glob = (fclaw2d_global_t*) user;
	fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
	const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

	fc2d_hps_patch patch;
	if (level == 0) {
		// Build root patch
		int Nx = clawpatch_opt->mx;
		int Ny = clawpatch_opt->my;
		double x_lower = fclaw_opt->ax;
		double x_upper = fclaw_opt->bx;
		double y_lower = fclaw_opt->ay;
		double y_upper = fclaw_opt->by;
		fc2d_hps_patchgrid root_grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);
		patch.grid = root_grid;
		patch.level = level;
		if (quadtree->child_indices[level][idx] == -1) {
			patch.is_leaf = true;
		}
		else {
			patch.is_leaf = false;
		}
		patch.N_cells_leaf = Nx;
		patch.N_patch_side = {1, 1, 1, 1};
		patch.has_coarsened = false;
		patch.user = (fclaw2d_global_t*) glob;
	}
	else {
		// Build branch patch
		// Get parent patch
		int pID = quadtree->parent_indices[level][idx];
		fc2d_hps_patch parent_patch = quadtree->data[pID];

		// Build child patch
		int Nx = parent_patch.grid.Nx;
		int Ny = parent_patch.grid.Ny;
		double x_midpoint = (parent_patch.grid.x_lower + parent_patch.grid.x_upper) / 2.0;
		double y_midpoint = (parent_patch.grid.y_lower + parent_patch.grid.y_upper) / 2.0;
		double x_lower, x_upper, y_lower, y_upper;
		if (idx % 4 == 0) {
			// Lower left
			x_lower = parent_patch.grid.x_lower;
			x_upper = x_midpoint;
			y_lower = parent_patch.grid.y_lower;
			y_upper = y_midpoint;
		}
		else if (idx % 4 == 1) {
			// Lower right
			x_lower = x_midpoint;
			x_upper = parent_patch.grid.x_upper;
			y_lower = parent_patch.grid.y_lower;
			y_upper = y_midpoint;
		}
		else if (idx % 4 == 2) {
			// Upper left
			x_lower = parent_patch.grid.x_lower;
			x_upper = x_midpoint;
			y_lower = y_midpoint;
			y_upper = parent_patch.grid.y_upper;
		}
		else if (idx % 4 == 3) {
			// Upper right
			x_lower = x_midpoint;
			x_upper = parent_patch.grid.x_upper;
			y_lower = y_midpoint;
			y_upper = parent_patch.grid.y_upper;
		}
		fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);
		patch.grid = grid;
		patch.level = level;
		if (quadtree->child_indices[level][idx] == -1) {
			patch.is_leaf = true;
		}
		else {
			patch.is_leaf = false;
		}
		patch.N_cells_leaf = Nx;
		patch.N_patch_side = {1, 1, 1, 1};
		patch.has_coarsened = false;
		patch.user = (fclaw2d_global_t*) glob;
	}

	return patch;
}

void visit_set_ID(fc2d_hps_patch& patch) {
    if (patch.is_leaf == true) {
        patch.ID = current_ID++;
    }
    else {
        patch.ID = -1;
    }
}

/**
 * Given a glob, setup the HPS method by building the quadtree with `fc2d_hps_patch`s
 */
void fc2d_hps_setup(struct fclaw2d_global* glob) {

    fclaw_global_essentialf("Begin HPS setup\n");

	// Get domain and p4est
	fclaw2d_domain_t* domain = glob->domain;
	p4est_wrap_t* p4est_wrap = (p4est_wrap_t*) domain->pp;
	p4est_t* p4est = p4est_wrap->p4est;

	// Get quadtree
	fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance(p4est, init_fn, glob);

    // Get HPS options
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);

    // Set patch IDs
    current_ID = 0;
    quadtree->traverse_postorder(visit_set_ID);

    // Set up DtN cache
    if (hps_opt->cache_T) {
        T_cache.resize(quadtree->global_indices.size());
    }

    fclaw_global_essentialf("End HPS setup\n");

    return;

}