#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps_patchsolver.hpp>
#include <fc2d_hps_poissonproblem.hpp>
#include <fc2d_hps_patchgrid.hpp>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>

TEST(PatchSolver, init) {

	fc2d_hps_patchsolver* solver;
	fc2d_hps_FISHPACK_solver FISHPACK_solver;
	// Other solvers go here...

	// Check all derived solvers for proper polymorphism
	// FISHPACK90
	solver = &FISHPACK_solver;
	EXPECT_EQ(solver->get_moniker(), "FISHPACK90");

}

TEST(FISHPACK, init) {

	fc2d_hps_FISHPACK_solver solver;
	EXPECT_EQ(solver.get_moniker(), "FISHPACK90");

}

TEST(FISHPACK, linear_solve) {

	// Set up grid
	int N_points_side = 16;
	int Nx = N_points_side;
	int Ny = N_points_side;
	double x_lower = -1;
	double x_upper = 1;
	double y_lower = -1;
	double y_upper = 1;
	int N_unknowns = Nx * Ny;
	int N_boundary_points = 2*Nx + 2*Ny;

	fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

	// Set up Poisson problem
	fc2d_hps_poisson_problem poisson(LINEAR, x_lower, x_upper, y_lower, y_upper);
	fc2d_hps_vector<double> f_data(N_unknowns);
	fc2d_hps_vector<double> u_data(N_unknowns);
	fc2d_hps_vector<double> g_west(Ny);
	fc2d_hps_vector<double> g_east(Ny);
	fc2d_hps_vector<double> g_south(Nx);
	fc2d_hps_vector<double> g_north(Nx);
	fc2d_hps_vector<double> g_data(N_boundary_points);

	// Fill Poisson data
	for (int i = 0; i < Nx; i++) {
		double x = grid.point(XDIM, i);
		
		for (int j = 0; j < Ny; j++) {
			double y = grid.point(YDIM, j);
			int running_index = i + j*Nx;

			f_data[running_index] = poisson.f(x, y);
			u_data[running_index] = poisson.u(x, y);
		}
	}

	// Fill boundary data
	for (int j = 0; j < Ny; j++) {
		g_west[j] = poisson.u(grid.x_lower, grid.point(YDIM, j));
		g_east[j] = poisson.u(grid.x_upper, grid.point(YDIM, j));
	}
	for (int i = 0; i < Nx; i++) {
		g_south[i] = poisson.u(grid.point(XDIM, i), grid.y_lower);
		g_north[i] = poisson.u(grid.point(XDIM, i), grid.y_upper);
	}
	g_data.intract(0*Nx, g_west);
	g_data.intract(1*Nx, g_east);
	g_data.intract(2*Nx, g_south);
	g_data.intract(3*Nx, g_north);

	// Exact solution


	// Print out data
	std::cout << "Solving Poisson Problem: " << poisson.ID << std::endl;
	printf("  Nx = %i\n", Nx);
	printf("  Ny = %i\n", Ny);
	printf("  x_lower = %8.4f    x_upper = %8.4f\n", grid.x_lower, grid.x_upper);
	printf("  y_lower = %8.4f    y_upper = %8.4f\n", grid.y_lower, grid.y_upper);
	printf("  Dirichlet Data:\n");
	for (int j = 0; j < Ny; j++) {
		printf("    g_west[%i]  = %8.4f    g_east[%i]  = %8.4f\n", j, g_west[j], j, g_east[j]);
	}
	std::cout << std::endl;
	for (int i = 0; i < Nx; i++) {
		printf("    g_south[%i] = %8.4f    g_north[%i] = %8.4f\n", i, g_south[i], i, g_north[i]);
	}
	printf("  Grid Data:\n");
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			printf("    x = %8.4f    y = %8.4f    rhs[%i, %i] = %8.4f    u[%i, %i] = %8.4f\n", grid.point(XDIM, i), grid.point(YDIM, j), i, j, f_data[j + i*Nx], i, j, u_data[j + i*Nx]);
		}
	}

	// Create solver and solve
	fc2d_hps_FISHPACK_solver solver;
	fc2d_hps_vector<double> u_FISHPACK = solver.solve(grid, g_data, f_data);

	// Print and check for correctness
	printf("  Solution Data:\n");
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			printf("    x = %8.4f    y = %8.4f    u_FIHSPACK[%i, %i] = %8.4f    u_EXACT[%i, %i] = %8.4f\n", grid.point(XDIM, i), grid.point(YDIM, j), i, j, u_FISHPACK[j + i*Nx], i, j, u_data[j + i*Nx]);
			EXPECT_NEAR(u_FISHPACK[j+i*Nx], u_data[j+i*Nx], 1e-14);
		}
	}

}

TEST(FISHPACK, solver_convergence) {

	std::vector<int> problems = {
		QUAD,
		POLY,
		TRIG	
	};

	for (auto& problem_ID : problems) {
		std::cout << "PROBLEM_ID = " << problem_ID << std::endl;

		int N_tests = 10;
		fc2d_hps_vector<double> errors(N_tests);
		fc2d_hps_vector<double> deltas(N_tests);
		for (int n = 0; n < N_tests; n++) {

			// Set up grid
			int N_points_side = pow(2, n+2);
			int Nx = N_points_side;
			int Ny = N_points_side;
			double x_lower = -1;
			double x_upper = 1;
			double y_lower = -1;
			double y_upper = 1;
			int N_unknowns = Nx * Ny;
			int N_boundary_points = 2*Nx + 2*Ny;

			fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

			// Set up Poisson problem
			fc2d_hps_poisson_problem poisson(problem_ID, x_lower, x_upper, y_lower, y_upper);
			fc2d_hps_vector<double> f_data(N_unknowns);
			fc2d_hps_vector<double> u_data(N_unknowns);
			fc2d_hps_vector<double> g_west(Ny);
			fc2d_hps_vector<double> g_east(Ny);
			fc2d_hps_vector<double> g_south(Nx);
			fc2d_hps_vector<double> g_north(Nx);
			fc2d_hps_vector<double> g_data(N_boundary_points);

			// Fill Poisson data
			for (int i = 0; i < Nx; i++) {
				double x = grid.point(XDIM, i);
				
				for (int j = 0; j < Ny; j++) {
					double y = grid.point(YDIM, j);
					// int running_index = i + j*Nx; // Transpose Order
					int running_index = j + i*Nx; // Patch Order

					f_data[running_index] = poisson.f(x, y);
					u_data[running_index] = poisson.u(x, y);
				}
			}

			// Fill boundary data
			for (int j = 0; j < Ny; j++) {
				g_west[j] = poisson.u(grid.x_lower, grid.point(YDIM, j));
				g_east[j] = poisson.u(grid.x_upper, grid.point(YDIM, j));
			}
			for (int i = 0; i < Nx; i++) {
				g_south[i] = poisson.u(grid.point(XDIM, i), grid.y_lower);
				g_north[i] = poisson.u(grid.point(XDIM, i), grid.y_upper);
			}
			g_data.intract(0*Nx, g_west);
			g_data.intract(1*Nx, g_east);
			g_data.intract(2*Nx, g_south);
			g_data.intract(3*Nx, g_north);

			// Create solver and solve
			fc2d_hps_FISHPACK_solver solver;
			fc2d_hps_vector<double> u_FISHPACK = solver.solve(grid, g_data, f_data);

			// Compute max error
			fc2d_hps_vector<double> u_error = u_FISHPACK - u_data;
			double max_diff = -1;
			for (int i = 0; i < N_unknowns; i++) {
				max_diff = fmax(max_diff, fabs(u_error[i]));
			}
			errors[n] = max_diff;
			deltas[n] = grid.dx;

			// Check for correctness
			printf("  N_unknowns = %8i    max_diff = %8.4e\n", N_unknowns, max_diff);
			
		}

		// Compute convergence order
		double order = 0;
		for (int n = 0; n < N_tests - 1; n++) {
			order = fmax(order, (log(errors[n]/errors[n+1])) / (log(deltas[n]/deltas[n+1])));
			printf("  order = %8.4e\n", order);
		}

		double expected_convergence_order = 1.8;
		EXPECT_GT(order, expected_convergence_order);
	}

}

TEST(FISHPACK, dtn_convergence) {

	std::vector<int> problems = {
		POLY,
		TRIG	
	};

	for (auto& problem_ID : problems) {
		std::cout << "PROBLEM_ID = " << problem_ID << std::endl;

		int N_tests = 10;
		fc2d_hps_vector<double> errors(N_tests);
		fc2d_hps_vector<double> deltas(N_tests);
		for (int n = 0; n < N_tests; n++) {

			// Set up grid
			int N_points_side = pow(2, n+2);
			int Nx = N_points_side;
			int Ny = N_points_side;
			double x_lower = -1;
			double x_upper = 1;
			double y_lower = -1;
			double y_upper = 1;
			int N_unknowns = Nx * Ny;
			int N_boundary_points = 2*Nx + 2*Ny;

			fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

			// Set up Poisson problem
			fc2d_hps_poisson_problem poisson(problem_ID, x_lower, x_upper, y_lower, y_upper);
			fc2d_hps_vector<double> f_data(N_unknowns);
			fc2d_hps_vector<double> g_west(Ny);
			fc2d_hps_vector<double> g_east(Ny);
			fc2d_hps_vector<double> g_south(Nx);
			fc2d_hps_vector<double> g_north(Nx);
			fc2d_hps_vector<double> g_data(N_boundary_points);
			fc2d_hps_vector<double> h_west(Ny);
			fc2d_hps_vector<double> h_east(Ny);
			fc2d_hps_vector<double> h_south(Nx);
			fc2d_hps_vector<double> h_north(Nx);
			fc2d_hps_vector<double> h_data(N_boundary_points);

			// Fill Poisson data
			for (int i = 0; i < Nx; i++) {
				double x = grid.point(XDIM, i);
				
				for (int j = 0; j < Ny; j++) {
					double y = grid.point(YDIM, j);
					// int running_index = i + j*Nx; // Transpose Order
					int running_index = j + i*Nx; // Patch Order

					f_data[running_index] = poisson.f(x, y);
				}
			}

			// Fill boundary data
			for (int j = 0; j < Ny; j++) {
				g_west[j] = poisson.u(grid.x_lower, grid.point(YDIM, j));
				g_east[j] = poisson.u(grid.x_upper, grid.point(YDIM, j));
				h_west[j] = poisson.dudx(grid.x_lower, grid.point(YDIM, j));
				h_east[j] = poisson.dudx(grid.x_upper, grid.point(YDIM, j));
			}
			for (int i = 0; i < Nx; i++) {
				g_south[i] = poisson.u(grid.point(XDIM, i), grid.y_lower);
				g_north[i] = poisson.u(grid.point(XDIM, i), grid.y_upper);
				h_south[i] = poisson.dudy(grid.point(XDIM, i), grid.y_lower);
				h_north[i] = poisson.dudy(grid.point(XDIM, i), grid.y_upper);
			}
			g_data.intract(0*Nx, g_west);
			g_data.intract(1*Nx, g_east);
			g_data.intract(2*Nx, g_south);
			g_data.intract(3*Nx, g_north);
			h_data.intract(0*Nx, h_west);
			h_data.intract(1*Nx, h_east);
			h_data.intract(2*Nx, h_south);
			h_data.intract(3*Nx, h_north);
			
			// Create solver and solve
			fc2d_hps_FISHPACK_solver solver;
			fc2d_hps_vector<double> h_FISHPACK = solver.dtn(grid, g_data, f_data);

			// Compute max error
			fc2d_hps_vector<double> h_error = h_FISHPACK - h_data;
			double max_diff = -1;
			for (int i = 0; i < 4*N_points_side; i++) {
				max_diff = fmax(max_diff, fabs(h_error[i]));
			}
			errors[n] = max_diff;
			deltas[n] = grid.dx;

			// Check for correctness
			printf("  N_unknowns = %8i    max_diff = %8.4e\n", N_unknowns, max_diff);
			
		}

		// Compute convergence order
		double order = 0;
		for (int n = 0; n < N_tests - 1; n++) {
			order = fmax(order, (log(errors[n]/errors[n+1])) / (log(deltas[n]/deltas[n+1])));
			printf("  order = %8.4e\n", order);
		}

		double expected_convergence_order = 0.8;
		EXPECT_GT(order, expected_convergence_order);
	}

}

TEST(FISHPACK, T_convergence) {

	std::vector<int> problems = {
		POLY,
		TRIG	
	};

	for (auto& problem_ID : problems) {
		std::cout << "PROBLEM_ID = " << problem_ID << std::endl;

		int N_tests = 7;
		fc2d_hps_vector<double> errors(N_tests);
		fc2d_hps_vector<double> deltas(N_tests);
		for (int n = 0; n < N_tests; n++) {

			// Set up grid
			int N_points_side = pow(2, n+2);
			int Nx = N_points_side;
			int Ny = N_points_side;
			double x_lower = -1;
			double x_upper = 1;
			double y_lower = -1;
			double y_upper = 1;
			int N_unknowns = Nx * Ny;
			int N_boundary_points = 2*Nx + 2*Ny;

			fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

			// Set up Poisson problem
			fc2d_hps_poisson_problem poisson(problem_ID, x_lower, x_upper, y_lower, y_upper);
			fc2d_hps_vector<double> f_data(N_unknowns);
			fc2d_hps_vector<double> g_west(Ny);
			fc2d_hps_vector<double> g_east(Ny);
			fc2d_hps_vector<double> g_south(Nx);
			fc2d_hps_vector<double> g_north(Nx);
			fc2d_hps_vector<double> g_data(N_boundary_points);
			fc2d_hps_vector<double> h_west(Ny);
			fc2d_hps_vector<double> h_east(Ny);
			fc2d_hps_vector<double> h_south(Nx);
			fc2d_hps_vector<double> h_north(Nx);
			fc2d_hps_vector<double> h_data(N_boundary_points);

			// Fill Poisson data
			for (int i = 0; i < Nx; i++) {
				double x = grid.point(XDIM, i);
				
				for (int j = 0; j < Ny; j++) {
					double y = grid.point(YDIM, j);
					// int running_index = i + j*Nx; // Transpose Order
					int running_index = j + i*Nx; // Patch Order

					f_data[running_index] = poisson.f(x, y);
				}
			}

			// Fill boundary data
			for (int j = 0; j < Ny; j++) {
				g_west[j] = poisson.u(grid.x_lower, grid.point(YDIM, j));
				g_east[j] = poisson.u(grid.x_upper, grid.point(YDIM, j));
				h_west[j] = poisson.dudx(grid.x_lower, grid.point(YDIM, j));
				h_east[j] = poisson.dudx(grid.x_upper, grid.point(YDIM, j));
			}
			for (int i = 0; i < Nx; i++) {
				g_south[i] = poisson.u(grid.point(XDIM, i), grid.y_lower);
				g_north[i] = poisson.u(grid.point(XDIM, i), grid.y_upper);
				h_south[i] = poisson.dudy(grid.point(XDIM, i), grid.y_lower);
				h_north[i] = poisson.dudy(grid.point(XDIM, i), grid.y_upper);
			}
			g_data.intract(0*Nx, g_west);
			g_data.intract(1*Nx, g_east);
			g_data.intract(2*Nx, g_south);
			g_data.intract(3*Nx, g_north);
			h_data.intract(0*Nx, h_west);
			h_data.intract(1*Nx, h_east);
			h_data.intract(2*Nx, h_south);
			h_data.intract(3*Nx, h_north);
			
			// Create solver and solve
			fc2d_hps_FISHPACK_solver solver;
			fc2d_hps_vector<double> h_FISHPACK = solver.dtn(grid, g_data, f_data);
			fc2d_hps_matrix<double> T_FISHPACK = solver.build_dtn(grid);

			// Build mapped Neumann data
			fc2d_hps_vector<double> g_zero(4*N_points_side, 0.0);
			fc2d_hps_vector<double> h_hat = solver.dtn(grid, g_zero, f_data);
			fc2d_hps_vector<double> h_mapped = T_FISHPACK * g_data;
			h_mapped = h_mapped + h_hat;

			// Compute max error
			fc2d_hps_vector<double> h_error = h_mapped - h_data;
			double max_diff = -1;
			for (int i = 0; i < 4*N_points_side; i++) {
				max_diff = fmax(max_diff, fabs(h_error[i]));
			}
			errors[n] = max_diff;
			deltas[n] = grid.dx;

			// Check for correctness
			printf("  N_unknowns = %8i    max_diff = %8.4e\n", N_unknowns, max_diff);
			
		}

		// Compute convergence order
		double order = 0;
		for (int n = 0; n < N_tests - 1; n++) {
			order = fmax(order, (log(errors[n]/errors[n+1])) / (log(deltas[n]/deltas[n+1])));
			printf("  order = %8.4e\n", order);
		}

		double expected_convergence_order = 0.8;
		EXPECT_GT(order, expected_convergence_order);
	}

}