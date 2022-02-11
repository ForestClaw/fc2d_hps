#include "Methods/fc2d_hps_setup.hpp"

patch_tree quadtree;
int current_ID;
std::vector<fc2d_hps_matrix<double>> T_cache;

void fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob) {

	// Build topmost patch
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
	const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
	int Nx = clawpatch_opt->mx;
	int Ny = clawpatch_opt->my;
	double x_lower = fclaw_opt->ax;
	double x_upper = fclaw_opt->bx;
	double y_lower = fclaw_opt->ay;
	double y_upper = fclaw_opt->by;
	fc2d_hps_patchgrid root_grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);
	fc2d_hps_patch root_patch;
	root_patch.grid = root_grid;
	root_patch.level = 0;
	root_patch.is_leaf = true;
	root_patch.N_cells_leaf = Nx;
	root_patch.N_patch_side = {1, 1, 1, 1};
	root_patch.user = (fclaw2d_global_t*) glob;

	// Create tree
    quadtree.root = new fc2d_hps_quadtree<fc2d_hps_patch>::fc2d_hps_quadnode(root_patch);
    quadtree.height = 1;

	// Build tree
	quadtree.build(build_from_p4est_callback_bigger, build_from_p4est_callback_init);
}

bool build_from_p4est_callback_bigger(fc2d_hps_patch& patch) {

	// Get access to p4est tree
	fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
	fclaw2d_domain_t* domain = glob->domain;
	p4est_wrap_t* p4est_wrap = (p4est_wrap_t*) domain->pp;
	p4est_t* p4est = p4est_wrap->p4est;
	p4est_tree_t* p4est_tree = p4est_tree_array_index(p4est->trees, 0); // TODO: What if this is part of a forest of trees?
	
	// Iterate over quadrants in p4est tree
	for (std::size_t i = 0; i < p4est->local_num_quadrants; i++) {

		// Get access to quadrant
		p4est_quadrant_t* p4est_quad = p4est_quadrant_array_index(&p4est_tree->quadrants, i);

		// Get levels of patch and quad
		int patch_level = patch.level;
		int quad_level = p4est_quad->level;

		if (quad_level > patch_level) {
			patch.is_leaf = false;
			return true;
		}
	}
	return false;

}

std::vector<fc2d_hps_patch> build_from_p4est_callback_init(fc2d_hps_patch& parent) {

	// Create children from parent patch
	//    Grid info
	int Nx = parent.grid.Nx;
	int Ny = parent.grid.Ny;
	double x_midpoint = (parent.grid.x_lower + parent.grid.x_upper) / 2;
	double y_midpoint = (parent.grid.y_lower + parent.grid.y_upper) / 2;

	//    Alpha (Lower Left)
	fc2d_hps_patchgrid alpha_grid(Nx, Ny, parent.grid.x_lower, x_midpoint, parent.grid.y_lower, y_midpoint);
	fc2d_hps_patch alpha_patch;
	alpha_patch.grid = alpha_grid;
	alpha_patch.level = parent.level + 1;
	alpha_patch.is_leaf = true;
	alpha_patch.N_cells_leaf = Nx;
	alpha_patch.N_patch_side = {1, 1, 1, 1};
	alpha_patch.user = (fclaw2d_global_t*) parent.user;

	//    Beta (Lower Right)
	fc2d_hps_patchgrid beta_grid(Nx, Ny, x_midpoint, parent.grid.x_upper, parent.grid.y_lower, y_midpoint);
	fc2d_hps_patch beta_patch;
	beta_patch.grid = beta_grid;
	beta_patch.level = parent.level + 1;
	beta_patch.is_leaf = true;
	beta_patch.N_cells_leaf = Nx;
	beta_patch.N_patch_side = {1, 1, 1, 1};
	beta_patch.user = (fclaw2d_global_t*) parent.user;

	//    Omega (Upper Left)
	fc2d_hps_patchgrid omega_grid(Nx, Ny, parent.grid.x_lower, x_midpoint, y_midpoint, parent.grid.y_upper);
	fc2d_hps_patch omega_patch;
	omega_patch.grid = omega_grid;
	omega_patch.level = parent.level + 1;
	omega_patch.is_leaf = true;
	omega_patch.N_cells_leaf = Nx;
	omega_patch.N_patch_side = {1, 1, 1, 1};
	omega_patch.user = (fclaw2d_global_t*) parent.user;

	//    Gamma (Upper Right)
	fc2d_hps_patchgrid gamma_grid(Nx, Ny, x_midpoint, parent.grid.x_upper, y_midpoint, parent.grid.y_upper);
	fc2d_hps_patch gamma_patch;
	gamma_patch.grid = gamma_grid;
	gamma_patch.level = parent.level + 1;
	gamma_patch.is_leaf = true;
	gamma_patch.N_cells_leaf = Nx;
	gamma_patch.N_patch_side = {1, 1, 1, 1};
	gamma_patch.user = (fclaw2d_global_t*) parent.user;

	std::vector<fc2d_hps_patch> children = {alpha_patch, beta_patch, omega_patch, gamma_patch};
	return children;
}

void visit_print2(fc2d_hps_patch& patch) {
    patch.print_info();
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
    fc2d_hps_create_quadtree_from_domain(glob);

    // Set patch IDs
    current_ID = 0;
    quadtree.traverse_inorder(visit_set_ID);

    // Set up DtN cache
    printf("setting up cache...\n");
    printf("quadtree.height = %i\n", quadtree.height);
    T_cache.reserve(quadtree.height);
    for (int i = 0; i < quadtree.height; i++) {
        T_cache[i] = fc2d_hps_matrix<double>(0,0);
    }
    // for (auto& T : T_cache) T = fc2d_hps_matrix<double>(0, 0);

    // T_cache_set.reserve(quadtree.height);

    printf("size of cache: %i\n", T_cache.size());

    fclaw_global_essentialf("End HPS setup\n");

    return;

}