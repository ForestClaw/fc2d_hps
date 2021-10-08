#include <iostream>
#include <fc2d_hps_patchsolver.hpp>
#include <fc2d_hps_poissonproblem.hpp>
#include <fc2d_hps_patchgrid.hpp>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>

int main() {
	std::cout << "Hello from test_bench" << std::endl;

	// Set up grid
	int N_points_side = 8;
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
			int running_index = j + i*Nx;

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
			// EXPECT_NEAR(u_FISHPACK[j+i*Nx], u_data[j+i*Nx], 1e-14);
		}
	}

	return 0;
}