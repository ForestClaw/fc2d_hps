#include <Methods/fc2d_hps_merge.hpp>

// Cache vector for DtN matrix
extern std::vector<fc2d_hps_matrix<double>> T_cache;

std::vector<int> fill_range(int start, int end) {
	std::vector<int> v(end - start);
	for (int i = 0; i < v.size(); i++) {
		v[i] = start + i;
	}
	return v;
}

fc2d_hps_matrix<double> build_L21(int n_rows, int n_cols) {
	fc2d_hps_matrix<double> L21(n_rows, n_cols, 0);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			if (i == 2*j) {
				L21(i,j) = 0.5;
				L21(i,j+1) = 0.5;
			}
		}
	}
	return L21;
}

fc2d_hps_matrix<double> build_L12(int n_rows, int n_cols) {
	fc2d_hps_matrix<double> L12(n_rows, n_cols, 0);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			if (i == 0 && j == 0) {
				L12(i,j) = 1.25;
				L12(i,j+1) = -0.25;
			}
			else if (i == n_rows-1 && j == n_cols-2) {
				L12(i,j) = 0.75;
				L12(i,j+1) = 0.25;
			}
			else if (j == 2*(i-1)) {
				L12(i,j) = 0.75;
				L12(i,j+1) = 0.25;
				L12(i+1,j) = 0.25;
				L12(i+1,j+1) = 0.75;
			}
		}
	}
	return L12;
}

fc2d_hps_matrix<double> merge_X(fc2d_hps_matrix<double>& T_33_alpha, fc2d_hps_matrix<double>& T_33_beta) {
	// X_inv = T_33_alpha - T_33_beta

	// std::cout << "[merge_X]  building X_tau" << std::endl;
	return T_33_alpha - T_33_beta;

}

fc2d_hps_matrix<double> merge_X_alpha_refined(fc2d_hps_matrix<double>& T_33_alpha, fc2d_hps_matrix<double>& T_33_beta, fc2d_hps_matrix<double>& L21, fc2d_hps_matrix<double>& L12) {
	fc2d_hps_matrix<double> temp = L21 * T_33_alpha;
	temp = temp * L12;
	return temp - T_33_beta;
}

fc2d_hps_matrix<double> merge_X_beta_refined(fc2d_hps_matrix<double>& T_33_alpha, fc2d_hps_matrix<double>& T_33_beta, fc2d_hps_matrix<double>& L21, fc2d_hps_matrix<double>& L12) {
	fc2d_hps_matrix<double> temp = L21 * T_33_beta;
	temp = temp * L12;
	return T_33_alpha - temp;
}

fc2d_hps_matrix<double> merge_S(fc2d_hps_matrix<double>& X_tau, fc2d_hps_matrix<double>& T_31_alpha, fc2d_hps_matrix<double>& T_32_beta) {
	// S = X_tau * [-T_31_alpha | T_32_beta] <-- Linear solve as X_tau is inverted

	// Form RHS
	// std::cout << "[merge_S]  building S_tau" << std::endl;
	fc2d_hps_matrix<double> RHS(T_31_alpha.rows, T_31_alpha.cols + T_32_beta.cols);
	T_31_alpha.negate();
	RHS.intract(0, 0, T_31_alpha);
	RHS.intract(0, T_31_alpha.cols, T_32_beta);
	T_31_alpha.negate();
	return solve(X_tau, RHS);

}

fc2d_hps_matrix<double> merge_S_alpha_refined(fc2d_hps_matrix<double>& X_tau, fc2d_hps_matrix<double>& T_31_alpha, fc2d_hps_matrix<double>& T_32_beta, fc2d_hps_matrix<double>& L21) {
	fc2d_hps_matrix<double> temp = L21 * T_31_alpha;
	temp.negate();
	fc2d_hps_matrix<double> temp2(temp.rows, temp.cols + T_32_beta.cols);
	return solve(X_tau, temp2);
}

fc2d_hps_matrix<double> merge_S_beta_refined(fc2d_hps_matrix<double>& X_tau, fc2d_hps_matrix<double>& T_31_alpha, fc2d_hps_matrix<double>& T_32_beta, fc2d_hps_matrix<double>& L21) {
	fc2d_hps_matrix<double> temp = L21 * T_32_beta;
	T_31_alpha.negate();
	fc2d_hps_matrix<double> temp2(T_31_alpha.rows, T_31_alpha.cols + temp.cols);
	return solve(X_tau, temp2);
}

fc2d_hps_matrix<double> merge_T(fc2d_hps_matrix<double>& S_tau, fc2d_hps_matrix<double>& T_11_alpha, fc2d_hps_matrix<double>& T_22_beta, fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta) {
	// T = [ [T_11_alpha, 0], [0, T_22_beta] ] + [ [T_13_alpha], [T_23_beta] ] * S_tau

	// Form left_term matrix
	// std::cout << "[merge_T]  building T_tau" << std::endl;
	fc2d_hps_matrix<double> left_term(T_11_alpha.rows + T_22_beta.rows, T_11_alpha.cols + T_22_beta.cols, 0.0);
	left_term.intract(0, 0, T_11_alpha);
	left_term.intract(T_11_alpha.rows, T_11_alpha.cols, T_22_beta);

	// Form right_term matrix
	fc2d_hps_matrix<double> H(T_13_alpha.rows + T_23_beta.rows, T_13_alpha.cols);
	H.intract(0, 0, T_13_alpha);
	H.intract(T_13_alpha.rows, 0, T_23_beta);
	fc2d_hps_matrix<double> right_term = H * S_tau;
	return left_term + right_term;

}

fc2d_hps_matrix<double> merge_T_alpha_refined(fc2d_hps_matrix<double>& S_tau, fc2d_hps_matrix<double>& T_11_alpha, fc2d_hps_matrix<double>& T_22_beta, fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_matrix<double>& L12) {
	fc2d_hps_matrix<double> left_term(T_11_alpha.rows + T_22_beta.rows, T_11_alpha.cols + T_22_beta.cols, 0.0);
	left_term.intract(0, 0, T_11_alpha);
	left_term.intract(T_11_alpha.rows, T_11_alpha.cols, T_22_beta);

	fc2d_hps_matrix<double> temp = T_13_alpha * L12;
	fc2d_hps_matrix<double> H(temp.rows + T_23_beta.rows, temp.cols);
	H.intract(0, 0, temp);
	H.intract(temp.rows, 0, T_23_beta);
	fc2d_hps_matrix<double> right_term = H * S_tau;

	return left_term + right_term;
}

fc2d_hps_matrix<double> merge_T_beta_refined(fc2d_hps_matrix<double>& S_tau, fc2d_hps_matrix<double>& T_11_alpha, fc2d_hps_matrix<double>& T_22_beta, fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_matrix<double>& L12) {
	fc2d_hps_matrix<double> left_term(T_11_alpha.rows + T_22_beta.rows, T_11_alpha.cols + T_22_beta.cols, 0.0);
	left_term.intract(0, 0, T_11_alpha);
	left_term.intract(T_11_alpha.rows, T_11_alpha.cols, T_22_beta);

	fc2d_hps_matrix<double> temp = T_23_beta * L12;
	fc2d_hps_matrix<double> H(T_13_alpha.rows + temp.rows, T_13_alpha.cols);
	H.intract(0, 0, T_13_alpha);
	H.intract(T_13_alpha.rows, 0, temp);
	fc2d_hps_matrix<double> right_term = H * S_tau;

	return left_term + right_term;
}

index_set_t make_index_sets_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

	int N_points_leaf_side = alpha.N_cells_leaf;
	int N_W_alpha = alpha.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_alpha = alpha.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_alpha = alpha.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_alpha = alpha.N_patch_side[NORTH] * N_points_leaf_side;
	int N_W_beta = beta.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_beta = beta.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_beta = beta.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_beta = beta.N_patch_side[NORTH] * N_points_leaf_side;

	std::vector<int> I_W_alpha = fill_range(0, N_W_alpha);
	std::vector<int> I_E_alpha = fill_range(N_W_alpha, N_W_alpha + N_E_alpha);
	std::vector<int> I_S_alpha = fill_range(N_W_alpha + N_E_alpha, N_W_alpha + N_E_alpha + N_S_alpha);
	std::vector<int> I_N_alpha = fill_range(N_W_alpha + N_E_alpha + N_S_alpha, N_W_alpha + N_E_alpha + N_S_alpha + N_N_alpha);
	std::vector<int> I_W_beta = fill_range(0, N_W_beta);
	std::vector<int> I_E_beta = fill_range(N_W_beta, N_W_beta + N_E_beta);
	std::vector<int> I_S_beta = fill_range(N_W_beta + N_E_beta, N_W_beta + N_E_beta + N_S_beta);
	std::vector<int> I_N_beta = fill_range(N_W_beta + N_E_beta + N_S_beta, N_W_beta + N_E_beta + N_S_beta + N_N_beta);
	std::vector<int> I_1(0);
	std::vector<int> I_2(0);
	std::vector<int> I_3_alpha(0);
	std::vector<int> I_3_beta(0);

	I_1.insert(I_1.end(), I_W_alpha.begin(), I_W_alpha.end());
	I_1.insert(I_1.end(), I_S_alpha.begin(), I_S_alpha.end());
	I_1.insert(I_1.end(), I_N_alpha.begin(), I_N_alpha.end());
	I_2.insert(I_2.end(), I_E_beta.begin(), I_E_beta.end());
	I_2.insert(I_2.end(), I_S_beta.begin(), I_S_beta.end());
	I_2.insert(I_2.end(), I_N_beta.begin(), I_N_beta.end());
	I_3_alpha.insert(I_3_alpha.end(), I_E_alpha.begin(), I_E_alpha.end());
	I_3_beta.insert(I_3_beta.end(), I_W_beta.begin(), I_W_beta.end());

	index_set_t index_set;
	index_set.I1 = I_1;
	index_set.I2 = I_2;
	index_set.I3_alpha = I_3_alpha;
	index_set.I3_beta = I_3_beta;

	return index_set;

}

fc2d_hps_patch merge_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {
	printf("horizontal merge: alpha ID = %i, beta ID = %i\n", alpha.ID, beta.ID);

	// Get options
	fclaw2d_global_t* glob = (fclaw2d_global_t*) alpha.user;
	fc2d_hps_options* hps_opt = fc2d_hps_get_options(glob);

	// Get number of points for each side for each patch
	int N_points_leaf_side = alpha.N_cells_leaf;
	int N_W_alpha = alpha.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_alpha = alpha.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_alpha = alpha.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_alpha = alpha.N_patch_side[NORTH] * N_points_leaf_side;
	int N_W_beta = beta.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_beta = beta.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_beta = beta.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_beta = beta.N_patch_side[NORTH] * N_points_leaf_side;

	// Build index vectors
	index_set_t index_sets = make_index_sets_horizontal(alpha, beta);

	// Extract blocks
	fc2d_hps_matrix<double> T_11_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I1);
	fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I3_alpha);
	fc2d_hps_matrix<double> T_31_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I1);
	fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I3_alpha);

	fc2d_hps_matrix<double> T_22_beta = beta.T.from_index_set(index_sets.I2, index_sets.I2);
	fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(index_sets.I2, index_sets.I3_beta);
	fc2d_hps_matrix<double> T_32_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I2);
	fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I3_beta);

	// Perform merge linear algebra
	// std::cout << "[merge_horizontal]  merging via linear algebra" << std::endl;
	fc2d_hps_matrix<double> X_tau;
	fc2d_hps_matrix<double> S_tau;
	fc2d_hps_matrix<double> T_tau;
	fc2d_hps_matrix<double> L21;
	fc2d_hps_matrix<double> L12;
	// Begin cases
	printf("alpha sides = [%i, %i, %i, %i]\n", alpha.N_patch_side[WEST], alpha.N_patch_side[EAST], alpha.N_patch_side[SOUTH], alpha.N_patch_side[NORTH]);
	printf("beta sides = [%i, %i, %i, %i]\n", beta.N_patch_side[WEST], beta.N_patch_side[EAST], beta.N_patch_side[SOUTH], beta.N_patch_side[NORTH]);
	if (alpha.N_patch_side[EAST] == beta.N_patch_side[WEST]) {
		// Uniform alpha and beta
		X_tau = merge_X(T_33_alpha, T_33_beta);
		S_tau = merge_S(X_tau, T_31_alpha, T_32_beta);
		T_tau = merge_T(S_tau, T_11_alpha, T_22_beta, T_13_alpha, T_23_beta);
	}
	else if (alpha.N_patch_side[EAST] == 2 * beta.N_patch_side[WEST]) {
		// Alpha refined
		printf("alpha refined\n");
		// printf("N_E_alpha = %i\n", N_E_alpha);
		// printf("N_W_beta = %i\n", N_W_beta);
		L21 = build_L21(N_W_beta, N_E_alpha);
		L12 = build_L12(N_E_alpha, N_W_beta);
		X_tau = merge_X_alpha_refined(T_33_alpha, T_33_beta, L21, L12);
		S_tau = merge_S_alpha_refined(X_tau, T_31_alpha, T_32_beta, L21);
		T_tau = merge_T_alpha_refined(S_tau, T_11_alpha, T_22_beta, T_13_alpha, T_23_beta, L12);

		// throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else if (2 * alpha.N_patch_side[EAST] == beta.N_patch_side[WEST]) {
		// Beta refined
		printf("beta refined\n");
		// printf("N_E_alpha = %i\n", N_E_alpha);
		// printf("N_W_beta = %i\n", N_W_beta);
		L12 = build_L12(N_W_beta, N_E_alpha);
		L21 = build_L21(N_E_alpha, N_W_beta);
		X_tau = merge_X_beta_refined(T_33_alpha, T_33_beta, L21, L12);
		S_tau = merge_S_beta_refined(X_tau, T_31_alpha, T_32_beta, L21);
		T_tau = merge_T_beta_refined(S_tau, T_11_alpha, T_22_beta, T_13_alpha, T_23_beta, L12);

		// throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else {
		throw std::invalid_argument("[fc2d_hps_merge::merge_horizontal] Size mismatch between alpha and beta.");
	}

	// Reorder matrices
	// std::cout << "[merge_horizontal]  reordering matrices" << std::endl;
	//    Reorder S
	std::vector<int> pi_nochange = {0};
	std::vector<int> pi_H = {0, 3, 1, 4, 2, 5};
	std::vector<int> R_S_tau = {N_E_alpha};
	std::vector<int> C_S_tau = {N_W_alpha, N_S_alpha, N_N_alpha, N_E_beta, N_S_beta, N_N_beta};
	printf("S.size = [%i, %i]\n", S_tau.rows, S_tau.cols);
	S_tau = S_tau.block_permute(pi_nochange, pi_H, R_S_tau, C_S_tau);

	//    Reorder T
	std::vector<int> R_T_tau = {N_W_alpha, N_S_alpha, N_N_alpha, N_E_beta, N_S_beta, N_N_beta};
	std::vector<int> C_T_tau = {N_W_alpha, N_S_alpha, N_N_alpha, N_E_beta, N_S_beta, N_N_beta};
	printf("T.size = [%i, %i]\n", T_tau.rows, T_tau.cols);
	T_tau = T_tau.block_permute(pi_H, pi_H, R_T_tau, C_T_tau);

	// Create new merged patch
	// std::cout << "[merge_horizontal]  creating merged patch" << std::endl;
	fc2d_hps_patchgrid merged_grid(alpha.grid.Nx + beta.grid.Nx, alpha.grid.Ny, alpha.grid.x_lower, beta.grid.x_upper, alpha.grid.y_lower, alpha.grid.y_upper);
	fc2d_hps_patch merged;
	merged.grid = merged_grid;
	merged.level = alpha.level;
	merged.is_leaf = false; 
	merged.N_patch_side = {
		alpha.N_patch_side[WEST],
		beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH] + beta.N_patch_side[SOUTH],
		alpha.N_patch_side[NORTH] + beta.N_patch_side[NORTH]
	};
	merged.N_cells_leaf = alpha.N_cells_leaf;
	merged.X = X_tau;
	merged.S = S_tau;
	merged.T = T_tau;

	// Set alpha to hold horizontal merge solution matrix
	alpha.T_prime = T_tau;
	alpha.S_prime = S_tau;

	// std::cout << "[merge_horizontal]  returning..." << std::endl;
	return merged;
	
}

index_set_t make_index_sets_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

	int N_points_leaf_side = alpha.N_cells_leaf;
	int N_W_alpha = alpha.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_alpha = alpha.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_alpha = alpha.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_alpha = alpha.N_patch_side[NORTH] * N_points_leaf_side;
	int N_W_beta = beta.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_beta = beta.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_beta = beta.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_beta = beta.N_patch_side[NORTH] * N_points_leaf_side;

	std::vector<int> I_W_alpha = fill_range(0, N_W_alpha);
	std::vector<int> I_E_alpha = fill_range(N_W_alpha, N_W_alpha + N_E_alpha);
	std::vector<int> I_S_alpha = fill_range(N_W_alpha + N_E_alpha, N_W_alpha + N_E_alpha + N_S_alpha);
	std::vector<int> I_N_alpha = fill_range(N_W_alpha + N_E_alpha + N_S_alpha, N_W_alpha + N_E_alpha + N_S_alpha + N_N_alpha);
	std::vector<int> I_W_beta = fill_range(0, N_W_beta);
	std::vector<int> I_E_beta = fill_range(N_W_beta, N_W_beta + N_E_beta);
	std::vector<int> I_S_beta = fill_range(N_W_beta + N_E_beta, N_W_beta + N_E_beta + N_S_beta);
	std::vector<int> I_N_beta = fill_range(N_W_beta + N_E_beta + N_S_beta, N_W_beta + N_E_beta + N_S_beta + N_N_beta);
	std::vector<int> I_1(0);
	std::vector<int> I_2(0);
	std::vector<int> I_3_alpha(0);
	std::vector<int> I_3_beta(0);

	I_1.insert(I_1.end(), I_W_alpha.begin(), I_W_alpha.end());
	I_1.insert(I_1.end(), I_E_alpha.begin(), I_E_alpha.end());
	I_1.insert(I_1.end(), I_S_alpha.begin(), I_S_alpha.end());
	I_2.insert(I_2.end(), I_W_beta.begin(), I_W_beta.end());
	I_2.insert(I_2.end(), I_E_beta.begin(), I_E_beta.end());
	I_2.insert(I_2.end(), I_N_beta.begin(), I_N_beta.end());
	I_3_alpha.insert(I_3_alpha.end(), I_N_alpha.begin(), I_N_alpha.end());
	I_3_beta.insert(I_3_beta.end(), I_S_beta.begin(), I_S_beta.end());

	index_set_t index_set;
	index_set.I1 = I_1;
	index_set.I2 = I_2;
	index_set.I3_alpha = I_3_alpha;
	index_set.I3_beta = I_3_beta;

	return index_set;

}

fc2d_hps_patch merge_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {
	printf("vertical merge: alpha ID = %i, beta ID = %i\n", alpha.ID, beta.ID);

	// Build index vectors
	int N_points_leaf_side = alpha.N_cells_leaf;
	int N_W_alpha = alpha.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_alpha = alpha.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_alpha = alpha.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_alpha = alpha.N_patch_side[NORTH] * N_points_leaf_side;
	int N_W_beta = beta.N_patch_side[WEST] * N_points_leaf_side;
	int N_E_beta = beta.N_patch_side[EAST] * N_points_leaf_side;
	int N_S_beta = beta.N_patch_side[SOUTH] * N_points_leaf_side;
	int N_N_beta = beta.N_patch_side[NORTH] * N_points_leaf_side;

	// Build index vectors
	index_set_t index_sets = make_index_sets_vertical(alpha, beta);

	// Extract blocks
	fc2d_hps_matrix<double> T_11_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I1);
	fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I3_alpha);
	fc2d_hps_matrix<double> T_31_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I1);
	fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I3_alpha);

	fc2d_hps_matrix<double> T_22_beta = beta.T.from_index_set(index_sets.I2, index_sets.I2);
	fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(index_sets.I2, index_sets.I3_beta);
	fc2d_hps_matrix<double> T_32_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I2);
	fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I3_beta);

	// Perform merge linear algebra
	fc2d_hps_matrix<double> X_tau;
	fc2d_hps_matrix<double> S_tau;
	fc2d_hps_matrix<double> T_tau;
	//    Begin cases
	//    Uniform alpha and beta
	if (alpha.N_patch_side[NORTH] == beta.N_patch_side[SOUTH]) {
		X_tau = merge_X(T_33_alpha, T_33_beta);
		S_tau = merge_S(X_tau, T_31_alpha, T_32_beta);
		T_tau = merge_T(S_tau, T_11_alpha, T_22_beta, T_13_alpha, T_23_beta);
	}
	// @TODO: Put in other cases
	else {
		throw std::invalid_argument("[fc2d_hps_merge::merge_vertical] Size mismatch between alpha and beta.");
	}

	// Reorder to get back to WESN ordering
	// std::cout << "[merge_vertical]  reordering matrices" << std::endl;
	std::vector<int> pi_nochange = {0};
	std::vector<int> pi_V = {0, 3, 1, 4, 2, 5};
	std::vector<int> R_S_tau = {N_N_alpha};
	std::vector<int> C_S_tau = {N_W_alpha, N_E_alpha, N_S_alpha, N_W_beta, N_E_beta, N_N_beta};
	S_tau = S_tau.block_permute(pi_nochange, pi_V, R_S_tau, C_S_tau);

	std::vector<int> R_T_tau = {N_W_alpha, N_E_alpha, N_S_alpha, N_W_beta, N_E_beta, N_N_beta};
	std::vector<int> C_T_tau = {N_W_alpha, N_E_alpha, N_S_alpha, N_W_beta, N_E_beta, N_N_beta};
	T_tau = T_tau.block_permute(pi_V, pi_V, R_T_tau, C_T_tau);

	// Create new merged patch
	// std::cout << "[merge_vertical]  creating merged patch" << std::endl;
	fc2d_hps_patchgrid merged_grid(alpha.grid.Nx, alpha.grid.Ny + beta.grid.Ny, alpha.grid.x_lower, alpha.grid.x_upper, alpha.grid.y_lower, beta.grid.y_upper);
	fc2d_hps_patch merged;
	merged.grid = merged_grid;
	merged.level = alpha.level - 1;
	merged.is_leaf = false;
	merged.N_patch_side = {
		alpha.N_patch_side[WEST] + beta.N_patch_side[WEST],
		alpha.N_patch_side[EAST] + beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH],
		beta.N_patch_side[NORTH]
	};
	merged.N_cells_leaf = alpha.N_cells_leaf;
	merged.X = X_tau;
	merged.S = S_tau;
	merged.T = T_tau;
	
	// std::cout << "[merge_vertical]  returning..." << std::endl;
	return merged;

}

void merge_4to1(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

	// Error checks
	//    Check if children data are set
	if (child0.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child0.T.size()` is 0; it shouldn't be..."); }
	if (child1.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child1.T.size()` is 0; it shouldn't be..."); }
	if (child2.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child2.T.size()` is 0; it shouldn't be..."); }
	if (child3.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child3.T.size()` is 0; it shouldn't be..."); }

	// Get options
	fclaw2d_global_t* glob = (fclaw2d_global_t*) parent.user;
	fc2d_hps_options* hps_opt = fc2d_hps_get_options(glob);

	// Forward declarations
	fc2d_hps_patch alpha_prime;
	fc2d_hps_patch beta_prime;
	fc2d_hps_patch tau;

	// if (hps_opt->cache_T) {
	// 	if (T_cache[parent.level].rows == 0 && T_cache[parent.level].cols == 0) {
	// 		// T for parent not set; compute and put in cache
	// 		alpha_prime = merge_horizontal(child0, child1);
	// 		beta_prime = merge_horizontal(child2, child3);
	// 		tau = merge_vertical(alpha_prime, beta_prime);
	// 		T_cache[parent.level] = tau.T;
	// 	}
	// }

	// Horizontal merge
	// std::cout << "[merge_4to1]  begin horizontal merge 1" << std::endl;
	alpha_prime = merge_horizontal(child0, child1);

	// std::cout << "[merge_4to1]  begin horizontal merge 2" << std::endl;
	beta_prime = merge_horizontal(child2, child3);

	// Vertical merge
	// std::cout << "[merge_4to1]  begin vertical merge" << std::endl;
	tau = merge_vertical(alpha_prime, beta_prime);
	
	// Copy only necessary data to parent
	// std::cout << "[merge_4to1]  begin copy to parent" << std::endl;
	parent.level = tau.level;
	parent.is_leaf = tau.is_leaf;
	parent.N_patch_side = tau.N_patch_side;
	parent.N_cells_leaf = tau.N_cells_leaf;
	parent.grid = tau.grid;
	parent.T = tau.T;
	parent.S = tau.S;
	parent.X = tau.X;

	// std::cout << "[merge_4to1]  end 4-to-1 merge, returning..." << std::endl;

}