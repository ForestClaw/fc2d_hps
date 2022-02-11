#include "Methods/fc2d_hps_copy.hpp"

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

    // std::string T_filename = "T_level_" + std::to_string(patch.level) + ".mmio";
    // std::FILE* T_file = fopen(T_filename.c_str(), "w");
    // patch.T.write_to_mmio(T_file, DataType::real, "%16.8e");
    // fclose(T_file);
}

void fc2d_hps_clawpatch_data_move(fclaw2d_global* glob) {
    fclaw_global_essentialf("Begin move to ForestClaw data\n");
    // fclaw_global_essentialf("!TODO!\n");

    quadtree.traverse_inorder(visit_copy_data);
    // hps_patch_quadtree_to_vtk(quadtree, "hps_patches");
    fclaw_global_essentialf("End move to ForestClaw data\n");
}