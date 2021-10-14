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

	// Transpose RHS data for FORTRAN call
	fc2d_hps_vector<double> fT(N_side * N_side);
	for (int i = 0; i < grid.Nx; i++) {
		for (int j = 0; j < grid.Ny; j++) {
			fT[i + j*N_side] = rhs_data[j + i*N_side];
		}
	}

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
	hstcrtt_(&A, &B, &M, &MBDCND, BDA, BDB,
			&C, &D, &N, &NBDCND, BDC, BDD,
			&ELMBDA, F, &IDIMF, &PERTRB, &IERROR, W);
	// hstcrt_(&A, &B, &M, &MBDCND, BDA, BDB,
	// 		&C, &D, &N, &NBDCND, BDC, BDD,
	// 		&ELMBDA, F, &IDIMF, &PERTRB, &IERROR);
	if (IERROR != 0) {
		std::cerr << "[fc2d_hps_FISHPACK_solver::solve] WARNING: call to hstcrt_ returned non-zero error value: IERROR = " << IERROR << std::endl;
	}

	// Move FISHPACK solution into fc2d_hps_vector for output
	fc2d_hps_vector<double> solution(grid.Nx * grid.Ny);
	for (int i = 0; i < grid.Nx; i++) {
		for (int j = 0; j < grid.Ny; j++) {
			solution[j + i*N_side] = F[i + j*N_side];
		}
	}

	return solution; // return rhs_data;

}

fc2d_hps_vector<double> fc2d_hps_FISHPACK_solver::dtn(fc2d_hps_patchgrid grid, fc2d_hps_vector<double> dirichlet_data, fc2d_hps_vector<double> rhs_data) {

	throw std::logic_error("[fc2d_hps_FISHPACK_solver::dtn] PLACEHOLDER; NOT IMPLEMENTED");
#if 0
	// Unpack grid data
	int N_cells = grid.Nx;

	// Unpack Dirichlet data
	fc2d_hps_vector<double> g_west = dirichlet_data.extract(0*N_cells, N_cells);
	fc2d_hps_vector<double> g_east = dirichlet_data.extract(1*N_cells, N_cells);
	fc2d_hps_vector<double> g_south = dirichlet_data.extract(2*N_cells, N_cells);
	fc2d_hps_vector<double> g_north = dirichlet_data.extract(3*N_cells, N_cells);

	// Compute solution on interior nodes
	fc2d_hps_vector<double> u = this->solve(grid, dirichlet_data, rhs_data);

	// Get interior edge cell data and compute Neumann data
	//    Interior cell data
	fc2d_hps_vector<double> u_west = u.extract();
	fc2d_hps_vector<double> u_east = u.extract();
	fc2d_hps_vector<double> u_south = u.extract();
	fc2d_hps_vector<double> u_north = u.extract();

	//    Neumann data
	double dtn_x = 2.0 / grid.dx;
	double dtn_y = 2.0 / grid.dy;
	fc2d_hps_vector<double> h_west = (u_west - g_west) * (dtn_x);
	fc2d_hps_vector<double> h_east = (u_east - g_east) * (-dtn_x);
	fc2d_hps_vector<double> h_south = (u_south - g_south) * (dtn_x);
	fc2d_hps_vector<double> h_north = (u_north - g_north) * (-dtn_x);

	//    Column stack and return
	fc2d_hps_vector<double> neumann_data(4*N_cells);
	neumann_data.intract(0*N_cells, h_west);
	neumann_data.intract(1*N_cells, h_east);
	neumann_data.intract(2*N_cells, h_south);
	neumann_data.intract(3*N_cells, h_north);
	return neumann_data;
#endif
	return dirichlet_data;
	
}

fc2d_hps_matrix<double> fc2d_hps_FISHPACK_solver::build_dtn(fc2d_hps_patchgrid grid) {

	throw std::logic_error("[fc2d_hps_FISHPACK_solver::build_dtn] PLACEHOLDER; NOT IMPLEMENTED");

	return fc2d_hps_matrix<double>{};
}