#ifndef FC2D_HPS_PATCHGRID_HPP
#define FC2D_HPS_PATCHGRID_HPP

#define XDIM 0		// Flag for x-dimension
#define YDIM 1		// Flag for y-dimension

#include <iostream>

class fc2d_hps_patchgrid {
	
public:

	// Grid variables
	int Nx, Ny;
	double x_lower, x_upper, y_lower, y_upper, dx, dy;

	// Constructors
	fc2d_hps_patchgrid();
	fc2d_hps_patchgrid(int Nx, int Ny, double x_lower, double x_upper, double y_lower, double y_upper);
	
	// Grid point accessor
	double point(int DIM, std::size_t index);

};

#endif // FC2D_HPS_PATCHGRID_HPP