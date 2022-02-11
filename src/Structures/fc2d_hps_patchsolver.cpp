#include <Structures/fc2d_hps_patchsolver.hpp>

#define DTN_OPTIMIZE 1

fc2d_hps_FISHPACK_solver::fc2d_hps_FISHPACK_solver() {}

std::string fc2d_hps_FISHPACK_solver::get_moniker() { return "FISHPACK90"; }

fc2d_hps_vector<double> fc2d_hps_FISHPACK_solver::solve(fc2d_hps_patchgrid grid, fc2d_hps_vector<double> dirichlet_data, fc2d_hps_vector<double> rhs_data) {

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

	std::size_t N = grid.Nx;
	std::size_t M = 4*N;
	fc2d_hps_matrix<double> T(M, M);
	fc2d_hps_vector<double> e_hat_j(M, 0.0);
	fc2d_hps_vector<double> f_zero(N*N, 0.0);
	fc2d_hps_vector<double> col_j(M);

#if DTN_OPTIMIZE
	// Iterate through first side of grid to form T
	// Compute first column of block T
	for (int j = 0; j < N; j++) {
		e_hat_j[j] = 1.0;
		col_j = this->dtn(grid, e_hat_j, f_zero);
		T.intract_column(j, col_j);
		e_hat_j[j] = 0.0;
	}

	// Extract blocks of T
	fc2d_hps_matrix<double> T_WW = T.extract(0*N, 0*N, N, N);
	fc2d_hps_matrix<double> T_EW = T.extract(1*N, 0*N, N, N);
	fc2d_hps_matrix<double> T_SW = T.extract(2*N, 0*N, N, N);
	fc2d_hps_matrix<double> T_NW = T.extract(3*N, 0*N, N, N);

	// Define other blocks in terms of first block column
	// T_WE = -T_EW
	// T_EE = -T_WW
	// T_SE = -T_NW^T
	// T_NE = Reversed columns from T_NW
	// 
	// T_WS = T_SW
	// T_ES = T_NW
	// T_SS = T_WW
	// T_NS = T_EW
	//
	// T_WN = -T_NW
	// T_EN = T_NE
	// T_SN = -T_EW
	// T_NN = -T_WW
	fc2d_hps_matrix<double> T_WE(N, N);
	fc2d_hps_matrix<double> T_EE(N, N);
	fc2d_hps_matrix<double> T_SE(N, N);
	fc2d_hps_matrix<double> T_NE(N, N);

	fc2d_hps_matrix<double> T_WS(N, N);
	fc2d_hps_matrix<double> T_ES(N, N);
	fc2d_hps_matrix<double> T_SS(N, N);
	fc2d_hps_matrix<double> T_NS(N, N);

	fc2d_hps_matrix<double> T_WN(N, N);
	fc2d_hps_matrix<double> T_EN(N, N);
	fc2d_hps_matrix<double> T_SN(N, N);
	fc2d_hps_matrix<double> T_NN(N, N);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			T_WE(i,j) = -T_EW(i,j);
			T_EE(i,j) = -T_WW(i,j);
			T_SE(i,j) = -T_NW(j,i);
			T_NE(i,j) = T_NW((N-1) - i, j);

			T_WS(i,j) = T_SW(i,j);
			T_ES(i,j) = T_NW(i,j);
			T_SS(i,j) = T_WW(i,j);
			T_NS(i,j) = T_EW(i,j);

			T_WN(i,j) = -T_NW(j,i);
			T_EN(i,j) = T_NE(i,j);
			T_SN(i,j) = -T_EW(i,j);
			T_NN(i,j) = -T_WW(i,j);
		}
	}

	// Intract blocks into T
	T.intract(0*N, 1*N, T_WE);
	T.intract(1*N, 1*N, T_EE);
	T.intract(2*N, 1*N, T_SE);
	T.intract(3*N, 1*N, T_NE);
	
	T.intract(0*N, 2*N, T_WS);
	T.intract(1*N, 2*N, T_ES);
	T.intract(2*N, 2*N, T_SS);
	T.intract(3*N, 2*N, T_NS);

	T.intract(0*N, 3*N, T_WN);
	T.intract(1*N, 3*N, T_EN);
	T.intract(2*N, 3*N, T_SN);
	T.intract(3*N, 3*N, T_NN);
#else
	// Iterate through all points on boundary to form T
	for (int j = 0; j < M; j++) {
		e_hat_j[j] = 1.0;
		col_j = this->dtn(grid, e_hat_j, f_zero);
		T.intract_column(j, col_j);
		e_hat_j[j] = 0.0;
	}
#endif
	return T;
}