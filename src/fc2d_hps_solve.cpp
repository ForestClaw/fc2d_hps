/*
Copyright (c) 2019-2021 Carsten Burstedde, Donna Calhoun, Damyn Chipman
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "fc2d_hps_solve.hpp"

fc2d_hps_quadtree<fc2d_hps_patch> fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob) {

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
	fc2d_hps_quadtree<fc2d_hps_patch> tree(root_patch);

	// Build tree
	tree.build(build_from_p4est_callback_bigger, build_from_p4est_callback_init);
	return tree;

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

/**
 * Given a glob, setup the HPS method by building the quadtree with `fc2d_hps_patch`s
 */
void fc2d_hps_setup(struct fclaw2d_global* glob) {

    fclaw_global_essentialf("Begin HPS setup\n");

    // Build quadtree and store in glob's user
    fc2d_hps_quadtree<fc2d_hps_patch> quadtree = fc2d_hps_create_quadtree_from_domain(glob);
    glob->user = (fc2d_hps_quadtree<fc2d_hps_patch>*) &quadtree;
    // printf("quadtree.root: %p\n", quadtree.root);
    // printf("glob: %p\n", glob);
    // printf("quadtree: %p\n", &quadtree);
    // printf("glob->user: %p\n", glob->user);

    fclaw_global_essentialf("End HPS setup\n");

    return;

}

/**
 * Callback function for a 4-to-1 merge
 */
static
void cb_merge(fclaw2d_domain_t *domain, fclaw2d_patch_t *fine_patches, int blockno, int fine0_patchno, void *user) {

    std::cout << "[fc2d_hps_solve.cpp::cb_merge]  In callback merge" << std::endl;
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
    // fc2d_hps_patch& tau = (fc2d_hps_patch*) glob->user;

    // Create four fc2d_hps_patchs from fclaw2d_patch_t*
    std::vector<fc2d_hps_patchgrid> grids(4); // @TODO: Redo to not have to build grids
    std::vector<fc2d_hps_patch> patches(4);
    for (int i = 0; i < 4; i++) {
        // Populate grid info
        int mbc;
        fclaw2d_clawpatch_grid_data(g->glob, &(fine_patches[i]), &(grids[i].Nx), &(grids[i].Ny), &mbc, &(grids[i].x_lower), &(grids[i].y_lower), &(grids[i].dx), &(grids[i].dy));
        grids[i].x_upper = grids[i].x_lower + grids[i].dx*grids[i].Nx;
        grids[i].y_upper = grids[i].y_lower + grids[i].dy*grids[i].Ny;

        // Create patches
        patches[i].grid = grids[i];
        patches[i].ID = blockno + i;
        // @TODO: Get current level
        patches[i].is_leaf = true;
        patches[i].N_patch_side[WEST] = 1;
        patches[i].N_patch_side[EAST] = 1;
        patches[i].N_patch_side[SOUTH] = 1;
        patches[i].N_patch_side[NORTH] = 1;
    }

    // Build DtNs for all patches
    fc2d_hps_FISHPACK_solver FISHPACK_solver;
    fc2d_hps_matrix<double> T = FISHPACK_solver.build_dtn(patches[0].grid);
    for (int i = 0; i < 4; i++) {
        patches[i].T = T;
    }

    // // Merge 4-to-1
    // merge_4to1(tau, patches[0], patches[1], patches[2], patches[3]);
    
}

void visit_print2(fc2d_hps_patch& patch) {
    patch.print_info();
}

void visit_leaves_build_dtn(fc2d_hps_patch& patch) {
    if (patch.is_leaf) {
        // printf("=================================================\n");
        fc2d_hps_FISHPACK_solver FISHPACK_solver;
        patch.T = FISHPACK_solver.build_dtn(patch.grid);
        // patch.print_info();
    }
    return;
}

void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    merge_4to1(tau, alpha, beta, gamma, omega);
}

void fc2d_hps_build(fclaw2d_global_t *glob)
{
    fclaw_global_essentialf("Begin HPS build\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch> quadtree = *(fc2d_hps_quadtree<fc2d_hps_patch>*) glob->user;

    // Build DtN matrix on leaf patches
    quadtree.traverse(visit_leaves_build_dtn);

    // Merge DtN's up tree and build solution
    quadtree.merge(visit_merge);

    fclaw_global_essentialf("End HPS build\n");
}

void visit_upwards(fc2d_hps_patch& patch) {
    
}

void fc2d_hps_upwards(fclaw2d_global_t* glob) {
    fclaw_global_essentialf("Begin HPS upwards pass\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch> quadtree = *(fc2d_hps_quadtree<fc2d_hps_patch>*) glob->user;

    // quadtree.traverse(visit_upwards);

    fclaw_global_essentialf("!TODO!\n");
    fclaw_global_essentialf("End HPS upwards pass\n");
}

void fc2d_hps_solve(fclaw2d_global_t* glob) {
    
    fclaw_global_essentialf("Begin HPS solve\n");

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    fclaw2d_domain_t* domain = glob->domain;
    
    /* Apply non-homogeneous boundary conditions to any patches on the boundary */
    if (fclaw_opt->maxlevel == 0) {
        fc2d_hps_physical_bc(glob);

        fclaw2d_patch_t* patch = &(domain->blocks->patches[0]);

        fc2d_hps_patchgrid grid;
        int mbc;
        fclaw2d_clawpatch_grid_data(glob, patch, &(grid.Nx), &(grid.Ny), &mbc, &(grid.x_lower), &(grid.y_lower), &(grid.dx), &(grid.dy));
        grid.x_upper = grid.x_lower + grid.dx*grid.Nx;
        grid.y_upper = grid.y_lower + grid.dy*grid.Ny;

        fc2d_hps_vector<double> g(4*grid.Nx, 0); // Pass in a bunch of zeros

        fc2d_hps_vector<double> f(grid.Nx * grid.Ny);
        int mfields;
        double* rhs;
        fclaw2d_clawpatch_rhs_data(glob, patch, &rhs, &mfields);
        for (int i = 0; i < (grid.Nx*grid.Ny); i++) {
            f[i] = rhs[i];

            printf("rhs[%i] = %16.8e\n", i, rhs[i]);
        }

        fc2d_hps_FISHPACK_solver solver;
        fc2d_hps_vector<double> u = solver.solve(grid, g, f);

        double* q;
        int meqn;
        fclaw2d_clawpatch_soln_data(glob, patch, &q, &meqn);
        printf("%p", q);

        int x_pts = grid.Nx + 2*mbc;
        int y_pts = grid.Ny + 2*mbc;
        for (int i = 0; i < x_pts; i++) {
            for (int j = 0; j < y_pts; j++) {
                if (i > 1 && i < x_pts-1 && j > 1 && j < y_pts-1) {
                    q[j + i*grid.Nx] = u[j + i*grid.Nx];
                    // rhs[j + i*grid.Nx] = u[j + i*grid.Nx];
                }
                else {
                    q[j + i*grid.Nx] = 0;
                    // rhs[j + i*grid.Nx] = 0;
                }
                printf("q[%i, %i] = %16.8e\n", i, j, q[j + i*grid.Nx]);
            }
        }
        for (int i = 0; i < (grid.Nx*grid.Ny); i++) {
            // q[i] = u[i];
        }
    }

    fclaw_global_essentialf("End HPS solve\n");
}
