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

/**
 * Callback function for a 4-to-1 merge
 */
static
void cb_merge(fclaw2d_domain_t *domain, fclaw2d_patch_t *fine_patches, int blockno, int fine0_patchno, void *user) {

    std::cout << "[fc2d_hps_solve.cpp::cb_merge]  In callback merge" << std::endl;
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
    // fc2d_hps_patch& tau = (fc2d_hps_patch*) glob->user;

    // Create four fc2d_hps_patchs from fclaw2d_patch_t*
    std::vector<fc2d_hps_patchgrid> grids(4); // @TOSO: Redo to not have to build grids
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

void fc2d_hps_build(fclaw2d_global_t *glob)
{
    std::cout << "[fc2d_hps_solve.cpp::fc2d_hps_build]  In build routine" << std::endl;
    // Get domain
    fclaw2d_domain* domain = glob->domain;

    /* Apply non-homogeneous boundary conditions to any patches on the boundary */
    // fc2d_hps_physical_bc(glob);

    /* Check that we only have one level */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    FCLAW_ASSERT(fclaw_opt->maxlevel <= 1);

    // Get options
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);

    // Call 
    // fc2d_hps_patch patch_tau_test;
    fclaw2d_global_iterate_families(glob, cb_merge, NULL); // @TODO: Change NULL to be a pointer to the parent patch
}

void fc2d_hps_solve(fclaw2d_global_t* glob) {
    std::cout << "[fc2d_hps_solve.cpp::fc2d_hps_solve]  In solve routine" << std::endl;

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
}
