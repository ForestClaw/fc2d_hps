#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <HPS/fc2d_hps.hpp>
#include <Structures/fc2d_hps_vector.hpp>
#include <Structures/fc2d_hps_matrix.hpp>
#include <Structures/fc2d_hps_patchsolver.hpp>
#include <Methods/fc2d_hps_merge.hpp>

TEST(Merge, merge_4to1_convergence) {

    std::cout << "Beginning unittest merge_4to1..." << std::endl;
    std::vector<int> Ns = {4, 8, 16, 32, 64, 128, 256};

    for (auto& N : Ns) {
        // Create 4 leaf patches
        //    Create four leaf grids
        std::cout << "  Building grids..." << std::endl;
        // int N = 16;
        fc2d_hps_patchgrid grid_alpha(N, N, -1.0, 0.0, -1.0, 0.0);
        fc2d_hps_patchgrid grid_beta(N, N, 0.0, 1.0, -1.0, 0.0);
        fc2d_hps_patchgrid grid_gamma(N, N, -1.0, 0.0, 0.0, 1.0);
        fc2d_hps_patchgrid grid_omega(N, N, 0.0, 1.0, 0.0, 1.0);

        //    Create four leaf patches
        std::cout << "  Building patches..." << std::endl;
        fc2d_hps_patch patch_alpha(grid_alpha, 0, 1, true);
        fc2d_hps_patch patch_beta(grid_beta, 1, 1, true);
        fc2d_hps_patch patch_gamma(grid_gamma, 2, 1, true);
        fc2d_hps_patch patch_omega(grid_omega, 3, 1, true);

        patch_alpha.N_cells_leaf = N;
        patch_beta.N_cells_leaf = N;
        patch_gamma.N_cells_leaf = N;
        patch_omega.N_cells_leaf = N;

        //    Build T on all patches
        std::cout << "  Building DtN's..." << std::endl;
        fc2d_hps_FISHPACK_solver FISHPACK_solver;
        patch_alpha.T = FISHPACK_solver.build_dtn(patch_alpha.grid);
        patch_beta.T = FISHPACK_solver.build_dtn(patch_beta.grid);
        patch_gamma.T = FISHPACK_solver.build_dtn(patch_gamma.grid);
        patch_omega.T = FISHPACK_solver.build_dtn(patch_omega.grid);

        //    Set h as zero (not testing this piece)
        patch_alpha.h = fc2d_hps_vector<double>(2*grid_alpha.Nx + 2*grid_alpha.Ny, 0);
        patch_beta.h = fc2d_hps_vector<double>(2*grid_beta.Nx + 2*grid_beta.Ny, 0);
        patch_gamma.h = fc2d_hps_vector<double>(2*grid_gamma.Nx + 2*grid_gamma.Ny, 0);
        patch_omega.h = fc2d_hps_vector<double>(2*grid_omega.Nx + 2*grid_omega.Ny, 0);

        // Compute parent patch
        std::cout << "  Performing 4-to-1 merge..." << std::endl;
        fc2d_hps_patch patch_tau_test;
        merge_4to1(patch_tau_test, patch_alpha, patch_beta, patch_gamma, patch_omega);
        
        // Create expected parent patch (refined single patch)
        //    Create parent grid
        fc2d_hps_patchgrid grid_tau(2*N, 2*N, -1.0, 1.0, -1.0, 1.0);

        //    Create parent patch
        fc2d_hps_patch patch_tau_expected(grid_tau, 4, 0, true);
        patch_tau_expected.N_cells_leaf = 2*N;

        //    Build T on parent patch
        patch_tau_expected.T = FISHPACK_solver.build_dtn(patch_tau_expected.grid);

        // Compare data
        EXPECT_EQ(patch_tau_expected.grid.Nx, patch_tau_test.grid.Nx);
        EXPECT_EQ(patch_tau_expected.grid.Ny, patch_tau_test.grid.Ny);
        EXPECT_EQ(patch_tau_expected.grid.x_lower, patch_tau_test.grid.x_lower);
        EXPECT_EQ(patch_tau_expected.grid.x_upper, patch_tau_test.grid.x_upper);
        EXPECT_EQ(patch_tau_expected.grid.y_lower, patch_tau_test.grid.y_lower);
        EXPECT_EQ(patch_tau_expected.grid.y_upper, patch_tau_test.grid.y_upper);
        EXPECT_EQ(patch_tau_expected.level, patch_tau_test.level);
        EXPECT_EQ(patch_tau_expected.T.rows, patch_tau_test.T.rows);
        EXPECT_EQ(patch_tau_expected.T.cols, patch_tau_test.T.cols);

        double error = 0;
        for (int i = 0; i < patch_tau_expected.T.rows; i++) {
            for (int j = 0; j < patch_tau_expected.T.cols; j++) {
                double diff = fabs(patch_tau_expected.T(i,j) - patch_tau_test.T(i,j));
                error = fmax(error, diff);
                // printf("i = %3i, j = %3i, E: %8.4f, T: %8.4f, diff = %8.4e, max_error = %8.4e\n", i, j, patch_tau_expected.T(i,j), patch_tau_test.T(i,j), diff, error);
            }
            std::cout << std::endl;
        }
        EXPECT_LT(error, 1e-12);
    }
}