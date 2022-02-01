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

// Global declarations
patch_tree quadtree; // use static...?
int current_ID;

/**================================================================================================
 * HPS Setup Routines
 *===============================================================================================*/
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

    // Build quadtree and store in glob's user
    // fc2d_hps_quadtree<fc2d_hps_patch>* quadtree_ptr = new fc2d_hps_quadtree<fc2d_hps_patch>;
    // fc2d_hps_quadtree<fc2d_hps_patch> quadtree = fc2d_hps_create_quadtree_from_domain(glob, interface);
    // fc2d_hps_quadtree<fc2d_hps_patch> fc2d_hps_interface::tree = fc2d_hps_create_quadtree_from_domain(glob);
    // fc2d_hps_interface::tree.traverse_inorder(visit_print2);
    // interface.tree = fc2d_hps_create_quadtree_from_domain(glob);
    // interface.tree.traverse_inorder(visit_print2);

    // fc2d_hps_quadtree<fc2d_hps_patch>* quadtree_ptr = &quadtree;
    // glob->user = quadtree_ptr;

    // printf("quadtree_ptr: %p\n", quadtree_ptr);
    // printf("quadtree_ptr->root: %p\n", quadtree_ptr->root);
    // printf("quadtree_ptr->root->children[0]: %p\n", quadtree_ptr->root->children[0]);
    // glob->user = (fc2d_hps_quadtree<fc2d_hps_patch>*) &quadtree;
    // printf("quadtree.root: %p\n", quadtree.root);
    // printf("glob: %p\n", glob);
    // printf("quadtree: %p\n", &quadtree);
    // printf("glob->user: %p\n", glob->user);

    fclaw_global_essentialf("End HPS setup\n");

    return;

}

/**================================================================================================
 * HPS Build Routines
 *===============================================================================================*/
void visit_leaves(fc2d_hps_patch& patch) {

    // TODO: Add level optimizations
    if (patch.is_leaf) {
        // Create patch solver
        fc2d_hps_FISHPACK_solver FISHPACK_solver;

        // Set DtN on leaf
        patch.T = FISHPACK_solver.build_dtn(patch.grid);

        // Set particular solution
        //    Get RHS data and set to f
        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fclaw2d_domain_t* domain = glob->domain;
        fclaw2d_patch_t* fc_patch = &(domain->blocks->patches[patch.ID]);
        int mfields;
        double* rhs;
        fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
        patch.f = fc2d_hps_vector<double>(rhs, rhs + patch.grid.Nx*patch.grid.Ny);

        //    Compute and set particular solution
        fc2d_hps_vector<double> g_zero(2*patch.grid.Nx + 2*patch.grid.Ny, 0);
        patch.w = FISHPACK_solver.solve(patch.grid, g_zero, patch.f);

        // Set Neumann data for particular solution
        patch.h = FISHPACK_solver.dtn(patch.grid, g_zero, patch.f);

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

    // Build DtN matrix on leaf patches
    quadtree.traverse_inorder(visit_leaves);

    // Merge DtN's up tree and build solution
    quadtree.merge(visit_merge);

    fclaw_global_essentialf("End HPS build\n");
}

/**================================================================================================
 * HPS Upwards Routines
 *===============================================================================================*/
void visit_upwards(fc2d_hps_patch& patch) {
    patch.print_info();
}

void fc2d_hps_upwards(fclaw2d_global_t* glob) {
    fclaw_global_essentialf("Begin HPS upwards pass\n");

    // Traverse inorder to set non-homogeneous RHS and BC data into solution vector
    // quadtree.traverse_inorder(visit_upwards);
    fclaw_global_essentialf("End HPS upwards pass\n");
}

/**================================================================================================
 * HPS Solve Routines
 *===============================================================================================*/
void visit_split(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega) {
    if (tau.is_leaf == false) {
        split_1to4(tau, alpha, beta, gamma, omega);
    }
}

void visit_patchsolver(fc2d_hps_patch& patch) {
    if (patch.is_leaf == true) {
        // Get patch solver
        // TODO: Put patch solver in glob...?
        fc2d_hps_FISHPACK_solver patch_solver;
        patch.u = patch_solver.solve(patch.grid, patch.g, patch.f);
    }
}

double qexact(double x, double y) {
    double b = (2.0/3.0) * M_PI;
    return sin(b*x) * sinh(b*y);
}

void fc2d_hps_solve(fclaw2d_global_t* glob) {
    
    fclaw_global_essentialf("Begin HPS solve\n");

    // Get quadtree
    // fc2d_hps_quadtree<fc2d_hps_patch> quadtree = *(fc2d_hps_quadtree<fc2d_hps_patch>*) glob->user;

    // Set top level Dirichlet boundary data (set to zero because ForestClaw handles it by moving it to RHS)
    // @TODO: Change this to set actual boundary conditions from glob
    int size_of_g = 2*quadtree.root->data.grid.Nx + 2*quadtree.root->data.grid.Ny;
    quadtree.root->data.g = fc2d_hps_vector<double>(size_of_g, 0);
    fc2d_hps_vector<double> g_west(quadtree.root->data.grid.Ny);
    fc2d_hps_vector<double> g_east(quadtree.root->data.grid.Ny);
    fc2d_hps_vector<double> g_south(quadtree.root->data.grid.Nx);
    fc2d_hps_vector<double> g_north(quadtree.root->data.grid.Nx);
    for (int j = 0; j < quadtree.root->data.grid.Ny; j++) {
        double y = quadtree.root->data.grid.point(YDIM, j);
        g_west[j] = qexact(quadtree.root->data.grid.x_lower, y);
        g_east[j] = qexact(quadtree.root->data.grid.x_upper, y);
    }
    for (int i = 0; i < quadtree.root->data.grid.Nx; i++) {
        double x = quadtree.root->data.grid.point(XDIM, i);
        g_south[i] = qexact(x, quadtree.root->data.grid.y_lower);
        g_north[i] = qexact(x, quadtree.root->data.grid.y_upper);
    }
    quadtree.root->data.g.intract(0*quadtree.root->data.grid.Nx, g_west);
    quadtree.root->data.g.intract(1*quadtree.root->data.grid.Nx, g_east);
    quadtree.root->data.g.intract(2*quadtree.root->data.grid.Nx, g_south);
    quadtree.root->data.g.intract(3*quadtree.root->data.grid.Nx, g_north);

    // Traverse tree from root and apply solution operator or patch solver
    quadtree.split(visit_split);

    // Iterate over leaf nodes and apply patch solver
    quadtree.traverse_inorder(visit_patchsolver);

    fclaw_global_essentialf("End HPS solve\n");
}

/**================================================================================================
 * HPS Interface Routines
 *===============================================================================================*/
void visit_copy_data(fc2d_hps_patch& patch) {
    if (patch.is_leaf == true) {

        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fclaw2d_domain_t* domain = glob->domain;
        fclaw2d_patch_t* fc_patch = &(domain->blocks->patches[patch.ID]);

        // Set pointer to solution data
        int mbc;
        int Nx, Ny;
        double x_lower, y_lower, dx, dy;
        double* q;
        int meqn, mfields;
        double* rhs;
        fclaw2d_clawpatch_grid_data(glob, fc_patch, &Nx, &Ny, &mbc, &x_lower, &y_lower, &dx, &dy);
        fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
        fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);

        // @TODO: Work on memcpy

        // Copy u into solution data
        int x_pts = patch.grid.Nx + 2*mbc;
        int y_pts = patch.grid.Ny + 2*mbc;
        for (int i = 0; i < x_pts; i++) {
            for (int j = 0; j < y_pts; j++) {
                int idx = j + i*y_pts;
                int idx_T = i + j*x_pts;
                if (i > mbc-1 && i < x_pts-mbc && j > mbc-1 && j < y_pts-mbc) {
                    q[idx_T] = patch.u[idx];
                    rhs[idx_T] = patch.u[idx];
                }
                else {
                    q[idx_T] = 0;
                    rhs[idx_T] = 0;
                }
                // printf("i = %i, j = %i, idx = %i, mbc = %i, x_pts = %i, y_pts = %i, q[%i] = %f\n", i, j, idx, mbc, x_pts, y_pts, idx, q[idx]);
            }
        }

        // Output patches
        // patch.to_vtk("patch_" + std::to_string(patch.ID), "", "wuf");
    }
}

void fc2d_hps_clawpatch_data_move(fclaw2d_global* glob) {
    fclaw_global_essentialf("Begin move to ForestClaw data\n");
    // fclaw_global_essentialf("!TODO!\n");

    quadtree.traverse_inorder(visit_copy_data);
    // hps_patch_quadtree_to_vtk(quadtree, "hps_patches");
    fclaw_global_essentialf("End move to ForestClaw data\n");
}