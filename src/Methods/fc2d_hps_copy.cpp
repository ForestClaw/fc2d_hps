#include "Methods/fc2d_hps_copy.hpp"
#include <stdio.h>

void visit_copy_data(fc2d_hps_patch& patch) {

    // TODO: Look at forestclaw iterator for copying

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
                    // q[idx_T] = 0;
                    // rhs[idx_T] = 0;
                }
            }
        }
    }

}

void visit_patch_to_vtk(fc2d_hps_patch& patch) {
    if (patch.is_leaf) {
        char filename_c[256];
        sprintf(filename_c, "patch%03i", patch.ID);
        std::string filename(filename_c);
        patch.to_vtk(filename, "Patch Data", "u");
    }
}

void fc2d_hps_clawpatch_data_move(fclaw2d_global* glob) {
    fclaw_global_essentialf("Begin move to ForestClaw data\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance();

    // quadtree->traverse_preorder(visit_patch_to_vtk);

    // Copy data to ForestClaw
    quadtree->traverse_preorder(visit_copy_data);
    
    fclaw_global_essentialf("End move to ForestClaw data\n");
}