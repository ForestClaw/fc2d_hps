#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_patchgrid.hpp>
#include <fc2d_hps_vector.hpp>

TEST(PatchGrid, init) {
	int Nx = 4;
	int Ny = 3;
	double x_lower = -1;
	double x_upper = 1;
	double y_lower = 0;
	double y_upper = 2;
	fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

	EXPECT_EQ(grid.Nx, Nx);
	EXPECT_EQ(grid.Ny, Ny);
	EXPECT_EQ(grid.x_lower, x_lower);
	EXPECT_EQ(grid.x_upper, x_upper);
	EXPECT_EQ(grid.y_lower, y_lower);
	EXPECT_EQ(grid.y_upper, y_upper);
}

TEST(PatchGrid, deltas) {
	int Nx = 4;
	int Ny = 3;
	double x_lower = -1;
	double x_upper = 1;
	double y_lower = 0;
	double y_upper = 2;
	double dx = (x_upper - x_lower)/Nx;
	double dy = (y_upper - y_lower)/Ny;
	fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

	EXPECT_EQ(grid.dx, dx);
	EXPECT_EQ(grid.dy, dy);
}

TEST(PatchGrid, points) {
	int Nx = 4;
	int Ny = 3;
	double x_lower = -1;
	double x_upper = 1;
	double y_lower = 0;
	double y_upper = 2;
	fc2d_hps_patchgrid grid(Nx, Ny, x_lower, x_upper, y_lower, y_upper);

	fc2d_hps_vector<double> x_points = {-3./4., -1./4., 1./4., 3./4.};
	fc2d_hps_vector<double> y_points = {1./3., 1., 5./3.};

	for (int i = 0; i < Nx; i++) {
		EXPECT_FLOAT_EQ(grid.point(XDIM, i), x_points[i]);
	}
	for (int j = 0; j < Ny; j++) {
		EXPECT_FLOAT_EQ(grid.point(YDIM, j), y_points[j]);
	}

	EXPECT_THROW(grid.point(XDIM, 4), std::out_of_range);
	EXPECT_THROW(grid.point(YDIM, 3), std::out_of_range);
	EXPECT_THROW(grid.point(XDIM, -2), std::out_of_range);
	EXPECT_THROW(grid.point(YDIM, -4), std::out_of_range);
	EXPECT_THROW(grid.point(2, 5), std::invalid_argument);
}