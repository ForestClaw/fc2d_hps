#include "gtest/gtest.h"
#include <vector>
#include <algorithm>
#include <string>
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>
#include <fc2d_hps_merge.hpp>
#include <fc2d_hps_split.hpp>
#include <fc2d_hps_patchsolver.hpp>
#include <fc2d_hps_poissonproblem.hpp>

TEST(Solve, solve_1to4_convergence) {

    // std::vector<int> Ns = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
    // std::vector<int> Ns = {4, 8, 16, 32, 64, 128};
    std::vector<int> Ns = {32};
    // std::vector<int> problem_IDs = {
    //     CONSTANT,
    //     LINEAR,
    //     LAPLACE1,
    //     QUAD,
    //     POLY,
    //     TRIG,
    //     GAUSSIAN
    // };
    std::vector<int> problem_IDs = {
        QUAD
    };
    std::vector<std::string> problem_names = {
        "CONSTANT : u = 1, f = 0",
        "LINEAR : u = x, f = 0",
        "LAPLACE : u = sin(b*x) * sinh(b*y), f = 0",
        "QUAD : u = x^2 + y^2 + 2*x*y, f = 4",
        "POLY : u = x^4 * y^2, f = 2*x^2 * (6*y^2 + x^2)",
        "TRIG : u = sin(2*PI*x) + sin(2*PI*y), f = -8*PI^2 * u",
        "GAUSSIAN : u = A*exp(-(x_part + y_part)), f = u * (stuff)"
    };

    bool debug_output = true;


    for (auto& problem_ID : problem_IDs) {

        printf("------------------------------------------ MAX ERRORS ------------------------------------------\n");
        printf(" PROBLEM ID : %i\n", problem_ID);
        printf(" %s\n", problem_names[problem_ID].c_str());
        printf(" ERROR TYPE     N             ALPHA         BETA          GAMMA         OMEGA         TOTAL\n");
        printf("------------  ------------  ------------  ------------  ------------  ------------  ------------\n");

        for (auto& N : Ns) {
            // Create Poisson problem
            // int problem_ID = TRIG;
            double x_lower = -1.0;
            double x_upper = 1.0;
            double y_lower = -1.0;
            double y_upper = 1.0;
            fc2d_hps_poisson_problem pde(problem_ID, x_lower, x_upper, y_lower, y_upper);

            // Create patch solver
            fc2d_hps_FISHPACK_solver FISHPACK_solver;
            
            // Create 4 leaf patches
            //    Create four leaf grids
            // int N = Ns[i];
            double x_halfway = (x_lower + x_upper) / 2.0;
            double y_halfway = (y_lower + y_upper) / 2.0;
            fc2d_hps_patchgrid grid_alpha(N, N, x_lower, x_halfway, y_lower, y_halfway);
            fc2d_hps_patchgrid grid_beta(N, N, x_halfway, x_upper, y_lower, y_halfway);
            fc2d_hps_patchgrid grid_gamma(N, N, x_lower, x_halfway, y_halfway, y_upper);
            fc2d_hps_patchgrid grid_omega(N, N, x_halfway, x_upper, y_halfway, y_upper);

            //    Create four leaf patches
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Building patches..." << std::endl;
            fc2d_hps_patch patch_alpha(grid_alpha, 0, 1, true);
            fc2d_hps_patch patch_beta(grid_beta, 1, 1, true);
            fc2d_hps_patch patch_gamma(grid_gamma, 2, 1, true);
            fc2d_hps_patch patch_omega(grid_omega, 3, 1, true);

            patch_alpha.N_cells_leaf = N;
            patch_beta.N_cells_leaf = N;
            patch_gamma.N_cells_leaf = N;
            patch_omega.N_cells_leaf = N;

            //    Build T on all patches
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Building DtN's..." << std::endl;
            patch_alpha.T = FISHPACK_solver.build_dtn(patch_alpha.grid);
            patch_beta.T = FISHPACK_solver.build_dtn(patch_beta.grid);
            patch_gamma.T = FISHPACK_solver.build_dtn(patch_gamma.grid);
            patch_omega.T = FISHPACK_solver.build_dtn(patch_omega.grid);

            //    Compute and set RHS data
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Computing RHS vectors `f`..." << std::endl;
            patch_alpha.f = fc2d_hps_vector<double>(grid_alpha.Nx*grid_alpha.Ny);
            for (int i = 0; i < grid_alpha.Nx; i++) {
                for (int j = 0; j < grid_alpha.Ny; j++) {
                    int idx = j + i*grid_alpha.Ny;
                    double x = grid_alpha.point(XDIM, i);
                    double y = grid_alpha.point(YDIM, j);
                    patch_alpha.f[idx] = pde.f(x, y);
                }
            }

            patch_beta.f = fc2d_hps_vector<double>(grid_beta.Nx*grid_beta.Ny);
            for (int i = 0; i < grid_beta.Nx; i++) {
                for (int j = 0; j < grid_beta.Ny; j++) {
                    int idx = j + i*grid_beta.Ny;
                    double x = grid_beta.point(XDIM, i);
                    double y = grid_beta.point(YDIM, j);
                    patch_beta.f[idx] = pde.f(x, y);
                }
            }

            patch_gamma.f = fc2d_hps_vector<double>(grid_gamma.Nx*grid_gamma.Ny);
            for (int i = 0; i < grid_gamma.Nx; i++) {
                for (int j = 0; j < grid_gamma.Ny; j++) {
                    int idx = j + i*grid_gamma.Ny;
                    double x = grid_gamma.point(XDIM, i);
                    double y = grid_gamma.point(YDIM, j);
                    patch_gamma.f[idx] = pde.f(x, y);
                }
            }

            patch_omega.f = fc2d_hps_vector<double>(grid_omega.Nx*grid_omega.Ny);
            for (int i = 0; i < grid_omega.Nx; i++) {
                for (int j = 0; j < grid_omega.Ny; j++) {
                    int idx = j + i*grid_omega.Ny;
                    double x = grid_omega.point(XDIM, i);
                    double y = grid_omega.point(YDIM, j);
                    patch_omega.f[idx] = pde.f(x, y);
                }
            }

            //    Compute particular solution on all patches
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Computing particular solution vectors `w`..." << std::endl;
            fc2d_hps_vector<double> g_zero(4*N, 0);
            patch_alpha.w = FISHPACK_solver.solve(grid_alpha, g_zero, patch_alpha.f);
            patch_beta.w = FISHPACK_solver.solve(grid_beta, g_zero, patch_beta.f);
            patch_gamma.w = FISHPACK_solver.solve(grid_gamma, g_zero, patch_gamma.f);
            patch_omega.w = FISHPACK_solver.solve(grid_omega, g_zero, patch_omega.f);

            //    Compute Neumann data for particular solution on all patches
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Computing Neumann data for particular solution `h`..." << std::endl;
            patch_alpha.h = FISHPACK_solver.dtn(grid_alpha, g_zero, patch_alpha.f);
            patch_beta.h = FISHPACK_solver.dtn(grid_beta, g_zero, patch_beta.f);
            patch_gamma.h = FISHPACK_solver.dtn(grid_gamma, g_zero, patch_gamma.f);
            patch_omega.h = FISHPACK_solver.dtn(grid_omega, g_zero, patch_omega.f);

            // patch_alpha.print_info();
            // patch_beta.print_info();
            // patch_gamma.print_info();
            // patch_omega.print_info();

            // Merge patches to one
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Merging patches..." << std::endl;
            fc2d_hps_patch patch_tau;
            merge_4to1(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);
            patch_tau.ID = -1;

            // Create Dirichlet data
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Creating parent Dirichlet data..." << std::endl;
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
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Beginning 1-to-4 solve..." << std::endl;
            patch_tau.g = g_tau;
            split_1to4(patch_tau, patch_alpha, patch_beta, patch_gamma, patch_omega);

            // Check that each leaf patch has accurate Dirichlet data
            double g_max_error = 0;
            std::vector<double> g_max_errors(4, 0);
            //    Alpha
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for alpha..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
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
                g_max_errors[0] = fmax(g_max_errors[0], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_errors[0]);
            }
            g_max_error = 0;

            //    Beta
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for beta..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
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
                g_max_errors[1] = fmax(g_max_errors[1], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_errors[1]);
            }
            g_max_error = 0;

            //    Gamma
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for gamma..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
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
                g_max_errors[2] = fmax(g_max_errors[2], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_errors[2]);
            }
            g_max_error = 0;

            //    Omega
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Dirichlet data for omega..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
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
                g_max_errors[3] = fmax(g_max_errors[3], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, g_exact, g_approx, diff, g_max_errors[3]);
            }


            // Compute Neumann data
            patch_alpha.v = patch_alpha.T * patch_alpha.g;
            patch_alpha.v = patch_alpha.v + patch_alpha.h;

            patch_beta.v = patch_beta.T * patch_beta.g;
            patch_beta.v = patch_beta.v + patch_beta.h;

            patch_gamma.v = patch_gamma.T * patch_gamma.g;
            patch_gamma.v = patch_gamma.v + patch_gamma.h;

            patch_omega.v = patch_omega.T * patch_omega.g;
            patch_omega.v = patch_omega.v + patch_omega.h;

            double v_max_error = 0;
            std::vector<double> v_max_errors(4, 0);
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Neumann data for alpha..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < 4*patch_alpha.grid.Nx; i++) {
                double x, y, v_exact;
                if (i >= 0*patch_alpha.grid.Nx && i < 1*patch_alpha.grid.Nx) {
                    // WEST
                    x = patch_alpha.grid.x_lower;
                    y = patch_alpha.grid.point(YDIM, i % patch_alpha.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 1*patch_alpha.grid.Nx && i < 2*patch_alpha.grid.Nx) {
                    // EAST
                    x = patch_alpha.grid.x_upper;
                    y = patch_alpha.grid.point(YDIM, i % patch_alpha.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 2*patch_alpha.grid.Nx && i < 3*patch_alpha.grid.Nx) {
                    // SOUTH
                    x = patch_alpha.grid.point(XDIM, i % patch_alpha.grid.Nx);
                    y = patch_alpha.grid.y_lower;
                    v_exact = pde.dudy(x, y);
                }
                else {
                    // NORTH
                    x = patch_alpha.grid.point(XDIM, i % patch_alpha.grid.Nx);
                    y = patch_alpha.grid.y_upper;
                    v_exact = pde.dudy(x, y);
                }

                double v_approx = patch_alpha.v[i];
                double diff = fabs(v_exact - v_approx);
                v_max_errors[0] = fmax(v_max_errors[0], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, v_exact, v_approx, diff, v_max_errors[0]);
            }
            v_max_error = 0;

            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Neumann data for beta..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < 4*patch_beta.grid.Nx; i++) {
                double x, y, v_exact;
                if (i >= 0*patch_beta.grid.Nx && i < 1*patch_beta.grid.Nx) {
                    // WEST
                    x = patch_beta.grid.x_lower;
                    y = patch_beta.grid.point(YDIM, i % patch_beta.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 1*patch_beta.grid.Nx && i < 2*patch_beta.grid.Nx) {
                    // EAST
                    x = patch_beta.grid.x_upper;
                    y = patch_beta.grid.point(YDIM, i % patch_beta.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 2*patch_beta.grid.Nx && i < 3*patch_beta.grid.Nx) {
                    // SOUTH
                    x = patch_beta.grid.point(XDIM, i % patch_beta.grid.Nx);
                    y = patch_beta.grid.y_lower;
                    v_exact = pde.dudy(x, y);
                }
                else {
                    // NORTH
                    x = patch_beta.grid.point(XDIM, i % patch_beta.grid.Nx);
                    y = patch_beta.grid.y_upper;
                    v_exact = pde.dudy(x, y);
                }

                double v_approx = patch_beta.v[i];
                double diff = fabs(v_exact - v_approx);
                v_max_errors[1] = fmax(v_max_errors[1], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, v_exact, v_approx, diff, v_max_errors[1]);
            }
            v_max_error = 0;

            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Neumann data for gamma..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < 4*patch_gamma.grid.Nx; i++) {
                double x, y, v_exact;
                if (i >= 0*patch_gamma.grid.Nx && i < 1*patch_gamma.grid.Nx) {
                    // WEST
                    x = patch_gamma.grid.x_lower;
                    y = patch_gamma.grid.point(YDIM, i % patch_gamma.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 1*patch_gamma.grid.Nx && i < 2*patch_gamma.grid.Nx) {
                    // EAST
                    x = patch_gamma.grid.x_upper;
                    y = patch_gamma.grid.point(YDIM, i % patch_gamma.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 2*patch_gamma.grid.Nx && i < 3*patch_gamma.grid.Nx) {
                    // SOUTH
                    x = patch_gamma.grid.point(XDIM, i % patch_gamma.grid.Nx);
                    y = patch_gamma.grid.y_lower;
                    v_exact = pde.dudy(x, y);
                }
                else {
                    // NORTH
                    x = patch_gamma.grid.point(XDIM, i % patch_gamma.grid.Nx);
                    y = patch_gamma.grid.y_upper;
                    v_exact = pde.dudy(x, y);
                }

                double v_approx = patch_gamma.v[i];
                double diff = fabs(v_exact - v_approx);
                v_max_errors[2] = fmax(v_max_errors[2], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, v_exact, v_approx, diff, v_max_errors[2]);
            }
            v_max_error = 0;

            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking Neumann data for omega..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < 4*patch_omega.grid.Nx; i++) {
                double x, y, v_exact;
                if (i >= 0*patch_omega.grid.Nx && i < 1*patch_omega.grid.Nx) {
                    // WEST
                    x = patch_omega.grid.x_lower;
                    y = patch_omega.grid.point(YDIM, i % patch_omega.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 1*patch_omega.grid.Nx && i < 2*patch_omega.grid.Nx) {
                    // EAST
                    x = patch_omega.grid.x_upper;
                    y = patch_omega.grid.point(YDIM, i % patch_omega.grid.Nx);
                    v_exact = pde.dudx(x, y);
                }
                else if (i >= 2*patch_omega.grid.Nx && i < 3*patch_omega.grid.Nx) {
                    // SOUTH
                    x = patch_omega.grid.point(XDIM, i % patch_omega.grid.Nx);
                    y = patch_omega.grid.y_lower;
                    v_exact = pde.dudy(x, y);
                }
                else {
                    // NORTH
                    x = patch_omega.grid.point(XDIM, i % patch_omega.grid.Nx);
                    y = patch_omega.grid.y_upper;
                    v_exact = pde.dudy(x, y);
                }

                double v_approx = patch_omega.v[i];
                double diff = fabs(v_exact - v_approx);
                v_max_errors[3] = fmax(v_max_errors[3], diff);

                if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, v_exact, v_approx, diff, v_max_errors[3]);
            }
            v_max_error = 0;
            
            // Use patchsolver to solve on each patch
            patch_alpha.u = FISHPACK_solver.solve(patch_alpha.grid, patch_alpha.g, patch_alpha.f);
            patch_beta.u = FISHPACK_solver.solve(patch_beta.grid, patch_beta.g, patch_beta.f);
            patch_gamma.u = FISHPACK_solver.solve(patch_gamma.grid, patch_gamma.g, patch_gamma.f);
            patch_omega.u = FISHPACK_solver.solve(patch_omega.grid, patch_omega.g, patch_omega.f);

            // Add particular solution...?
            // patch_alpha.u = patch_alpha.u + patch_alpha.w;
            // patch_beta.u = patch_beta.u + patch_beta.w;
            // patch_gamma.u = patch_gamma.u + patch_gamma.w;
            // patch_omega.u = patch_omega.u + patch_omega.w;

            // Fill tau patch with exact solution for comparison
            patch_tau.u = fc2d_hps_vector<double>(patch_tau.grid.Nx*patch_tau.grid.Ny);
            for (int i = 0; i < patch_tau.grid.Nx; i++) {
                for (int j = 0; j < patch_tau.grid.Ny; j++) {
                    double x = patch_tau.grid.point(XDIM, i);
                    double y = patch_tau.grid.point(YDIM, j);
                    int idx = j + i*patch_tau.grid.Ny;
                    patch_tau.u[idx] = pde.u(x, y);
                }
            }

            // Compare numerical solution to exact
            double max_error = 0;
            std::vector<double> max_errors(4, 0);
            //    Alpha
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking alpha..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < patch_alpha.grid.Nx; i++) {
                for (int j = 0; j < patch_alpha.grid.Ny; j++) {
                    double x = patch_alpha.grid.point(XDIM, i);
                    double y = patch_alpha.grid.point(YDIM, j);
                    double exact = pde.u(x,y);
                    double approx = patch_alpha.u[j + i*patch_alpha.grid.Nx];
                    double diff = fabs(exact - approx);
                    max_errors[0] = fmax(max_errors[0], diff);

                    if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_errors[0]);
                }
            }

            //    Beta
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking beta..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < patch_beta.grid.Nx; i++) {
                for (int j = 0; j < patch_beta.grid.Ny; j++) {
                    double x = patch_beta.grid.point(XDIM, i);
                    double y = patch_beta.grid.point(YDIM, j);
                    double exact = pde.u(x,y);
                    double approx = patch_beta.u[j + i*patch_beta.grid.Nx];
                    double diff = fabs(exact - approx);
                    max_errors[1] = fmax(max_errors[1], diff);

                    if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_errors[1]);
                }
            }

            //    Gamma
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking gamma..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < patch_gamma.grid.Nx; i++) {
                for (int j = 0; j < patch_gamma.grid.Ny; j++) {
                    double x = patch_gamma.grid.point(XDIM, i);
                    double y = patch_gamma.grid.point(YDIM, j);
                    double exact = pde.u(x,y);
                    double approx = patch_gamma.u[j + i*patch_gamma.grid.Nx];
                    double diff = fabs(exact - approx);
                    max_errors[2] = fmax(max_errors[2], diff);

                    if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_errors[2]);
                }
            }

            //    Omega
            if (debug_output) std::cout << "[unittest_solve::solve_1to4_convergence]  Checking omega..." << std::endl;
            if (debug_output) printf("  X             Y             EXACT         APPROX        DIFF          MAX_ERROR\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------\n");
            for (int i = 0; i < patch_omega.grid.Nx; i++) {
                for (int j = 0; j < patch_omega.grid.Ny; j++) {
                    double x = patch_omega.grid.point(XDIM, i);
                    double y = patch_omega.grid.point(YDIM, j);
                    double exact = pde.u(x,y);
                    double approx = patch_omega.u[j + i*patch_omega.grid.Nx];
                    double diff = fabs(exact - approx);
                    max_errors[3] = fmax(max_errors[3], diff);

                    if (debug_output) printf("%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", x, y, exact, approx, diff, max_errors[3]);
                }
            }

            g_max_error = *max_element(g_max_errors.begin(), g_max_errors.end());
            v_max_error = *max_element(v_max_errors.begin(), v_max_errors.end());
            max_error = *max_element(max_errors.begin(), max_errors.end());

            // patch_alpha.print_info();
            // patch_beta.print_info();
            // patch_gamma.print_info();
            // patch_omega.print_info();
            // patch_tau.print_info();

            if (debug_output) printf("------------------------------------------ MAX ERRORS ------------------------------------------\n");
            if (debug_output) printf(" ERROR TYPE     N             ALPHA         BETA          GAMMA         OMEGA         TOTAL\n");
            if (debug_output) printf("------------  ------------  ------------  ------------  ------------  ------------  ------------\n");
            printf(" DIRICHLET  | %12i %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", N, g_max_errors[0], g_max_errors[1], g_max_errors[2], g_max_errors[3], g_max_error);
            printf(" NEUMANN    | %12i %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", N, v_max_errors[0], v_max_errors[1], v_max_errors[2], v_max_errors[3], v_max_error);
            printf(" SOLUTION   | %12i %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n", N, max_errors[0], max_errors[1], max_errors[2], max_errors[3], max_error);

            if (N <= 32) {
                patch_alpha.to_vtk("alpha", "alpha patch data for problem: " + problem_names[problem_ID], "fwu");
                patch_beta.to_vtk("beta", "beta patch data for problem: " + problem_names[problem_ID], "fwu");
                patch_gamma.to_vtk("gamma", "gamma patch data for problem: " + problem_names[problem_ID], "fwu");
                patch_omega.to_vtk("omega", "omega patch data for problem: " + problem_names[problem_ID], "fwu");
                patch_tau.to_vtk("tau", "tau patch data for problem: " + problem_names[problem_ID], "u");
            }
        }
    }

    // EXPECT_LT(max_error, 1e-2);
}