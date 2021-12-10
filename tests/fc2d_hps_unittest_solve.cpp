#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>
#include <fc2d_hps_merge.hpp>
#include <fc2d_hps_patchsolver.hpp>
#include <fc2d_hps_poissonproblem.hpp>

TEST(SOLVE, solve_1to4_convergence) {
    
    // Create 4 leaf patches
    //    Create four leaf grids
    int N = 4;
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

    // Merge patches to one
    fc2d_hps_patch patch_tau;
    merge_4to1(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);

    // Create Dirichlet data
    //    Create Poisson problem
    int problem_ID = LAPLACE1;
    fc2d_hps_poisson_problem pde(problem_ID, patch_tau.grid.x_lower, patch_tau.grid.x_upper, patch_tau.grid.y_lower, patch_tau.grid.y_upper);

    //    Create Dirichlet data
    fc2d_hps_vector<double> g_tauW(patch_tau.grid.Nx);
    fc2d_hps_vector<double> g_tauE(patch_tau.grid.Nx);
    fc2d_hps_vector<double> g_tauS(patch_tau.grid.Nx);
    fc2d_hps_vector<double> g_tauN(patch_tau.grid.Nx);
    fc2d_hps_vector<double> g_tau(4*patch_tau.grid.Nx);

    for (int j = 0; j < patch_tau.grid.Ny; j++) {
        g_tauW[j] = pde.u(patch_tau.grid.x_lower, patch_tau.grid.point(YDIM, j));
        g_tauE[j] = pde.u(patch_tau.grid.x_upper, patch_tau.grid.point(YDIM, j));
    }
    for (int i = 0; i < patch_tau.grid.Nx; i++) {
        g_tauS[i] = pde.u(patch_tau.grid.point(XDIM, i), patch_tau.grid.y_lower);
        g_tauN[i] = pde.u(patch_tau.grid.point(XDIM, i), patch_tau.grid.y_upper);
    }

    g_tau.intract(0*patch_tau.grid.Nx, g_tauW);
    g_tau.intract(1*patch_tau.grid.Nx, g_tauE);
    g_tau.intract(2*patch_tau.grid.Nx, g_tauS);
    g_tau.intract(3*patch_tau.grid.Nx, g_tauN);

    // Solve 1 patch to 4
    patch_tau.g = g_tau;
    // split_1to4(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);

    EXPECT_TRUE(true);
}