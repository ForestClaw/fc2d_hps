#include "fc2d_hps_patchsolver.hpp"

fc2d_hps_FISHPACK_solver::fc2d_hps_FISHPACK_solver() {}

std::string fc2d_hps_FISHPACK_solver::get_moniker() { return "FISHPACK90"; }

fc2d_hps_vector<double> fc2d_hps_FISHPACK_solver::solve(fc2d_hps_patchgrid grid, fc2d_hps_vector<double> dirichlet_data, fc2d_hps_vector<double> rhs_data) {
	
	// throw std::logic_error("[fc2d_hps_FISHPACK_solver::solve] PLACEHOLDER; NOT IMPLEMENTED");

	// Unpack Dirichlet Data
	int N_side = grid.Nx;
	fc2d_hps_vector<double> g_west = dirichlet_data.extract(0*N_side, N_side);
	fc2d_hps_vector<double> g_east = dirichlet_data.extract(1*N_side, N_side);
	fc2d_hps_vector<double> g_south = dirichlet_data.extract(2*N_side, N_side);
	fc2d_hps_vector<double> g_north = dirichlet_data.extract(3*N_side, N_side);

	// Transpose RHS for FORTRAN call
	// @NOTE: Does this need to be done prior to function call as a vector?

	// Setup FORTRAN call to FISHPACK
	double A = grid.x_lower;
	double B = grid.x_upper;
	int M = grid.Nx;
	int MBDCND = 1;
	double* BDA = g_west.data();
	double* BDB = g_east.data();
	double C = grid.y_lower;
	double D = grid.y_upper;
	int N = grid.Ny;
	int NBDCND = 1;
	double* BDC = g_south.data();
	double* BDD = g_north.data();
	double ELMBDA = 0; // @TODO: Implement or get lambda value
	double* F = rhs_data.data();
	int IDIMF = M;
	double PERTRB;
	int IERROR;
	int WSIZE = 13*M + 4*N + M*((int)log2(N));
	double* W = (double*) malloc(WSIZE*sizeof(double));

	// Make FORTRAN call to FISHPACK
	// std::cout << "Calling hstcrt..." << std::endl;
	hstcrtt_(&A, &B, &M, &MBDCND, BDA, BDB,
			&C, &D, &N, &NBDCND, BDC, BDD,
			&ELMBDA, F, &IDIMF, &PERTRB, &IERROR, W);
	// hstcrt_(&A, &B, &M, &MBDCND, BDA, BDB,
	// 		&C, &D, &N, &NBDCND, BDC, BDD,
	// 		&ELMBDA, F, &IDIMF, &PERTRB, &IERROR);
	// std::cout << "Done with hstcrt!" << std::endl;
	if (IERROR != 0) {
		std::cerr << "[fc2d_hps_FISHPACK_solver::solve] WARNING: call to hstcrt_ returned non-zero error value: IERROR = " << IERROR << std::endl;
	}

	// Move FISHPACK solution into fc2d_hps_vector for output
	// @TODO
	fc2d_hps_vector<double> solution(grid.Nx * grid.Ny);
	for (int i = 0; i < grid.Nx; i++) {
		for (int j = 0; j < grid.Ny; j++) {
			int running_index = j + i*grid.Nx;
			solution[running_index] = F[running_index];
		}
	}

	return solution; // return rhs_data;

}

fc2d_hps_vector<double> fc2d_hps_FISHPACK_solver::dtn(fc2d_hps_patchgrid grid, fc2d_hps_vector<double> dirichlet_data, fc2d_hps_vector<double> rhs_data) {

	throw std::logic_error("[fc2d_hps_FISHPACK_solver::dtn] PLACEHOLDER; NOT IMPLEMENTED");

	return dirichlet_data;
}

fc2d_hps_matrix<double> fc2d_hps_FISHPACK_solver::build_dtn(fc2d_hps_patchgrid grid) {

	throw std::logic_error("[fc2d_hps_FISHPACK_solver::build_dtn] PLACEHOLDER; NOT IMPLEMENTED");

	return fc2d_hps_matrix<double>{};
}