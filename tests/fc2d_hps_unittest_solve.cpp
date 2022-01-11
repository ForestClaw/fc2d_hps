#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>
#include <fc2d_hps_merge.hpp>
#include <fc2d_hps_split.hpp>
#include <fc2d_hps_patchsolver.hpp>
#include <fc2d_hps_poissonproblem.hpp>

TEST(Solve, solve_1to4_convergence) {
    
    // Create 4 leaf patches
    //    Create four leaf grids
    int N = 16;
    fc2d_hps_patchgrid grid_alpha(N, N, -1.0, 0.0, -1.0, 0.0);
    fc2d_hps_patchgrid grid_beta(N, N, 0.0, 1.0, -1.0, 0.0);
    fc2d_hps_patchgrid grid_gamma(N, N, -1.0, 0.0, 0.0, 1.0);
    fc2d_hps_patchgrid grid_omega(N, N, 0.0, 1.0, 0.0, 1.0);

    //    Create four leaf patches
    std::cout << "[unittest_solve::solve_1to4_convergence]  Building patches..." << std::endl;
    fc2d_hps_patch patch_alpha(grid_alpha, 0, 1, true);
    fc2d_hps_patch patch_beta(grid_beta, 1, 1, true);
    fc2d_hps_patch patch_gamma(grid_gamma, 2, 1, true);
    fc2d_hps_patch patch_omega(grid_omega, 3, 1, true);

    patch_alpha.N_cells_leaf = N;
    patch_beta.N_cells_leaf = N;
    patch_gamma.N_cells_leaf = N;
    patch_omega.N_cells_leaf = N;

    //    Build T on all patches
    std::cout << "[unittest_solve::solve_1to4_convergence]  Building DtN's..." << std::endl;
    fc2d_hps_FISHPACK_solver FISHPACK_solver;
    patch_alpha.T = FISHPACK_solver.build_dtn(patch_alpha.grid);
    patch_beta.T = FISHPACK_solver.build_dtn(patch_beta.grid);
    patch_gamma.T = FISHPACK_solver.build_dtn(patch_gamma.grid);
    patch_omega.T = FISHPACK_solver.build_dtn(patch_omega.grid);

    // Merge patches to one
    std::cout << "[unittest_solve::solve_1to4_convergence]  Merging patches..." << std::endl;
    fc2d_hps_patch patch_tau;
    merge_4to1(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);

    // Create Dirichlet data
    std::cout << "[unittest_solve::solve_1to4_convergence]  Creating parent Dirichlet data..." << std::endl;
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
    std::cout << "[unittest_solve::solve_1to4_convergence]  Beginning 1-to-4 solve..." << std::endl;
    patch_tau.g = g_tau;
    split_1to4(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);

    // Check that each leaf patch has accurate Dirichlet data
    double g_max_error = 0;
    //    Alpha
    std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for alpha..." << std::endl;
    printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < 4*patch_alpha.grid.Nx; i++) {
        double x, y;
        if (i >= 0*patch_alpha.grid.Nx && i < 1*patch_alpha.grid.Nx) {
            // WEST
            x = patch_alpha.grid.x_lower;
            y = patch_alpha.grid.point(YDIM, i % patch_alpha.grid.Nx);
        }
        else if (i >= 1*patch_alpha.grid.Nx && i < 2*patch_alpha.grid.Nx) {
            // EAST
            x = patch_alpha.grid.x_upper;
            y = patch_alpha.grid.point(YDIM, i % patch_alpha.grid.Nx);
        }
        else if (i >= 2*patch_alpha.grid.Nx && i < 3*patch_alpha.grid.Nx) {
            // SOUTH
            x = patch_alpha.grid.point(XDIM, i % patch_alpha.grid.Nx);
            y = patch_alpha.grid.y_lower;
        }
        else {
            // NORTH
            x = patch_alpha.grid.point(XDIM, i % patch_alpha.grid.Nx);
            y = patch_alpha.grid.y_upper;
        }

        double g_approx = patch_alpha.g[i];
        double g_exact = pde.u(x, y);
        double diff = fabs(g_exact - g_approx);
        g_max_error = fmax(g_max_error, diff);

        printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_error);
    }
    g_max_error = 0;

    //    Beta
    std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for beta..." << std::endl;
    printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < 4*patch_beta.grid.Nx; i++) {
        double x, y;
        if (i >= 0*patch_beta.grid.Nx && i < 1*patch_beta.grid.Nx) {
            // WEST
            x = patch_beta.grid.x_lower;
            y = patch_beta.grid.point(YDIM, i % patch_beta.grid.Nx);
        }
        else if (i >= 1*patch_beta.grid.Nx && i < 2*patch_beta.grid.Nx) {
            // EAST
            x = patch_beta.grid.x_upper;
            y = patch_beta.grid.point(YDIM, i % patch_beta.grid.Nx);
        }
        else if (i >= 2*patch_beta.grid.Nx && i < 3*patch_beta.grid.Nx) {
            // SOUTH
            x = patch_beta.grid.point(XDIM, i % patch_beta.grid.Nx);
            y = patch_beta.grid.y_lower;
        }
        else {
            // NORTH
            x = patch_beta.grid.point(XDIM, i % patch_beta.grid.Nx);
            y = patch_beta.grid.y_upper;
        }

        double g_approx = patch_beta.g[i];
        double g_exact = pde.u(x, y);
        double diff = fabs(g_exact - g_approx);
        g_max_error = fmax(g_max_error, diff);

        printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_error);
    }
    g_max_error = 0;

    //    Gamma
    std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for gamma..." << std::endl;
    printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < 4*patch_gamma.grid.Nx; i++) {
        double x, y;
        if (i >= 0*patch_gamma.grid.Nx && i < 1*patch_gamma.grid.Nx) {
            // WEST
            x = patch_gamma.grid.x_lower;
            y = patch_gamma.grid.point(YDIM, i % patch_gamma.grid.Nx);
        }
        else if (i >= 1*patch_gamma.grid.Nx && i < 2*patch_gamma.grid.Nx) {
            // EAST
            x = patch_gamma.grid.x_upper;
            y = patch_gamma.grid.point(YDIM, i % patch_gamma.grid.Nx);
        }
        else if (i >= 2*patch_gamma.grid.Nx && i < 3*patch_gamma.grid.Nx) {
            // SOUTH
            x = patch_gamma.grid.point(XDIM, i % patch_gamma.grid.Nx);
            y = patch_gamma.grid.y_lower;
        }
        else {
            // NORTH
            x = patch_gamma.grid.point(XDIM, i % patch_gamma.grid.Nx);
            y = patch_gamma.grid.y_upper;
        }

        double g_approx = patch_gamma.g[i];
        double g_exact = pde.u(x, y);
        double diff = fabs(g_exact - g_approx);
        g_max_error = fmax(g_max_error, diff);

        printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_error);
    }
    g_max_error = 0;

    //    Omega
    std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for omega..." << std::endl;
    printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < 4*patch_omega.grid.Nx; i++) {
        double x, y;
        if (i >= 0*patch_omega.grid.Nx && i < 1*patch_omega.grid.Nx) {
            // WEST
            x = patch_omega.grid.x_lower;
            y = patch_omega.grid.point(YDIM, i % patch_omega.grid.Nx);
        }
        else if (i >= 1*patch_omega.grid.Nx && i < 2*patch_omega.grid.Nx) {
            // EAST
            x = patch_omega.grid.x_upper;
            y = patch_omega.grid.point(YDIM, i % patch_omega.grid.Nx);
        }
        else if (i >= 2*patch_omega.grid.Nx && i < 3*patch_omega.grid.Nx) {
            // SOUTH
            x = patch_omega.grid.point(XDIM, i % patch_omega.grid.Nx);
            y = patch_omega.grid.y_lower;
        }
        else {
            // NORTH
            x = patch_omega.grid.point(XDIM, i % patch_omega.grid.Nx);
            y = patch_omega.grid.y_upper;
        }

        double g_approx = patch_omega.g[i];
        double g_exact = pde.u(x, y);
        double diff = fabs(g_exact - g_approx);
        g_max_error = fmax(g_max_error, diff);

        printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_error);
    }

    // Set leaf patch RHS data to zero
    patch_alpha.f = fc2d_hps_vector<double>(N*N, 0);
    patch_beta.f = fc2d_hps_vector<double>(N*N, 0);
    patch_gamma.f = fc2d_hps_vector<double>(N*N, 0);
    patch_omega.f = fc2d_hps_vector<double>(N*N, 0);

    // Use patchsolver to solve on each patch
    patch_alpha.u = FISHPACK_solver.solve(patch_alpha.grid, patch_alpha.g, patch_alpha.f);
    patch_beta.u = FISHPACK_solver.solve(patch_beta.grid, patch_beta.g, patch_beta.f);
    patch_gamma.u = FISHPACK_solver.solve(patch_gamma.grid, patch_gamma.g, patch_gamma.f);
    patch_omega.u = FISHPACK_solver.solve(patch_omega.grid, patch_omega.g, patch_omega.f);

    // Compare numerical solution to exact
    double max_error = 0;
    //    Alpha
    std::cout << "[unittest::solve_1to4_convergence]  Checking alpha..." << std::endl;
    // printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    // printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < patch_alpha.grid.Nx; i++) {
        for (int j = 0; j < patch_alpha.grid.Ny; j++) {
            double x = patch_alpha.grid.point(XDIM, i);
            double y = patch_alpha.grid.point(YDIM, j);
            double exact = pde.u(x,y);
            double approx = patch_alpha.u[j + i*patch_alpha.grid.Nx];
            double diff = fabs(exact - approx);
            max_error = fmax(max_error, diff);

            // printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_error);
        }
    }

    //    Beta
    std::cout << "[unittest::solve_1to4_convergence]  Checking beta..." << std::endl;
    // printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    // printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < patch_beta.grid.Nx; i++) {
        for (int j = 0; j < patch_beta.grid.Ny; j++) {
            double x = patch_beta.grid.point(XDIM, i);
            double y = patch_beta.grid.point(YDIM, j);
            double exact = pde.u(x,y);
            double approx = patch_beta.u[j + i*patch_beta.grid.Nx];
            double diff = fabs(exact - approx);
            max_error = fmax(max_error, diff);

            // printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_error);
        }
    }

    //    Gamma
    std::cout << "[unittest::solve_1to4_convergence]  Checking gamma..." << std::endl;
    // printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    // printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < patch_gamma.grid.Nx; i++) {
        for (int j = 0; j < patch_gamma.grid.Ny; j++) {
            double x = patch_gamma.grid.point(XDIM, i);
            double y = patch_gamma.grid.point(YDIM, j);
            double exact = pde.u(x,y);
            double approx = patch_gamma.u[j + i*patch_gamma.grid.Nx];
            double diff = fabs(exact - approx);
            max_error = fmax(max_error, diff);

            // printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_error);
        }
    }

    //    Omega
    std::cout << "[unittest::solve_1to4_convergence]  Checking omega..." << std::endl;
    // printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
    // printf("------------  ------------  ------------  ------------  ------------  ------------\n");
    for (int i = 0; i < patch_omega.grid.Nx; i++) {
        for (int j = 0; j < patch_omega.grid.Ny; j++) {
            double x = patch_omega.grid.point(XDIM, i);
            double y = patch_omega.grid.point(YDIM, j);
            double exact = pde.u(x,y);
            double approx = patch_omega.u[j + i*patch_omega.grid.Nx];
            double diff = fabs(exact - approx);
            max_error = fmax(max_error, diff);

            // printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_error);
        }
    }

    patch_alpha.print_info();
    patch_beta.print_info();
    patch_gamma.print_info();
    patch_omega.print_info();
    patch_tau.print_info();

    EXPECT_TRUE(false);
}