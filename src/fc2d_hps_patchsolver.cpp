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
	double* F = fT.data();
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
		// std::cerr << "[fc2d_hps_FISHPACK_solver::solve] WARNING: call to hstcrt_ returned non-zero error value: IERROR = " << IERROR << std::endl;
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

	// throw std::logic_error("[fc2d_hps_FISHPACK_solver::dtn] PLACEHOLDER; NOT IMPLEMENTED");

	// Unpack grid data
	int N_side = grid.Nx;

	// Unpack Dirichlet data
	fc2d_hps_vector<double> g_west = dirichlet_data.extract(0*N_side, N_side);
	fc2d_hps_vector<double> g_east = dirichlet_data.extract(1*N_side, N_side);
	fc2d_hps_vector<double> g_south = dirichlet_data.extract(2*N_side, N_side);
	fc2d_hps_vector<double> g_north = dirichlet_data.extract(3*N_side, N_side);

	// Compute solution on interior nodes
	fc2d_hps_vector<double> u = this->solve(grid, dirichlet_data, rhs_data);

	// Get interior edge cell data and compute Neumann data
	//    Interior cell data
	fc2d_hps_vector<double> u_west(N_side);
	fc2d_hps_vector<double> u_east(N_side);
	fc2d_hps_vector<double> u_south(N_side);
	fc2d_hps_vector<double> u_north(N_side);
	
	//    Fill interior cell data
	for (int j = 0; j < N_side; j++) {
		u_west[j] = u[j];
		u_east[j] = u[(N_side-1)*N_side + j];
	}
	for (int i = 0; i < N_side; i++) {
		u_south[i] = u[i*N_side];
		u_north[i] = u[(i+1)*N_side - 1];
	}
	// if (N_side <= 4) {
	// 	for (int i = 0; i < N_side; i++) {
	// 		for (int j = 0; j < N_side; j++) {
	// 			printf("u[%i, %i] = u[%i] = %8.4f    ", i, j, i*N_side + j, u[i*N_side + j]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	for (int i = 0; i < N_side; i++) {
	// 		printf("u_west[%i] = %8.4f    u_east[%i] = %8.4f    u_south[%i] = %8.4f    u_north[%i] = %8.4f\n",
	// 				i, u_west[i], i, u_east[i], i, u_south[i], i, u_north[i]);
	// 	}
	// }

	//    Neumann data
	double dtn_x = 2.0 / grid.dx;
	double dtn_y = 2.0 / grid.dy;
	fc2d_hps_vector<double> h_west(N_side);
	fc2d_hps_vector<double> h_east(N_side);
	fc2d_hps_vector<double> h_south(N_side);
	fc2d_hps_vector<double> h_north(N_side);
	for (int i = 0; i < N_side; i++) {
		h_west[i]  = (dtn_x)  * (u_west[i] - g_west[i]);
		h_east[i]  = (-dtn_x) * (u_east[i] - g_east[i]);
		h_south[i] = (dtn_y)  * (u_south[i] - g_south[i]);
		h_north[i] = (-dtn_y) * (u_north[i] - g_north[i]);
	}

	//    Column stack and return
	fc2d_hps_vector<double> neumann_data(4*N_side);
	neumann_data.intract(0*N_side, h_west);
	neumann_data.intract(1*N_side, h_east);
	neumann_data.intract(2*N_side, h_south);
	neumann_data.intract(3*N_side, h_north);
	return neumann_data;
	
}

fc2d_hps_matrix<double> fc2d_hps_FISHPACK_solver::build_dtn(fc2d_hps_patchgrid grid) {

	throw std::logic_error("[fc2d_hps_FISHPACK_solver::build_dtn] PLACEHOLDER; NOT IMPLEMENTED");

	return fc2d_hps_matrix<double>{};
}