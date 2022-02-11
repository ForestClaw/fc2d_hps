#include <Structures/fc2d_hps_patchgrid.hpp>

fc2d_hps_patchgrid::fc2d_hps_patchgrid() :
	Nx(0), Ny(0), x_lower(0), x_upper(0), y_lower(0), y_upper(0), dx(0), dy(0)
		{}

fc2d_hps_patchgrid::fc2d_hps_patchgrid(int Nx, int Ny, double x_lower, double x_upper, double y_lower, double y_upper) :
	Nx(Nx), Ny(Ny), x_lower(x_lower), x_upper(x_upper), y_lower(y_lower), y_upper(y_upper), dx((x_upper - x_lower)/Nx), dy((y_upper - y_lower)/Ny)
		{}

double fc2d_hps_patchgrid::point(int DIM, std::size_t index) {
	if (DIM == XDIM) {
		if (index >= Nx || index < 0) {
			throw std::out_of_range("[fc2d_hps_patchgrid::point] `index` is out of range");
		}
		return (x_lower + dx/2) + index*dx;
	}
	else if (DIM == YDIM) {
		if (index >= Ny || index < 0) {
			throw std::out_of_range("[fc2d_hps_patchgrid::point] `index` is out of range");
		}
		return (y_lower + dy/2) + index*dy;
	}
	else {
		throw std::invalid_argument("[fc2d_hps_patchgrid::point] `DIM` is not a correct index");
	}
}
	