#include <Methods/fc2d_hps_upwards.hpp>

fc2d_hps_vector<double> merge_w(fc2d_hps_matrix<double>& X_tau, fc2d_hps_vector<double>& h_3_alpha, fc2d_hps_vector<double>& h_3_beta) {
	
	// Compute w
	fc2d_hps_vector<double> flux = h_3_beta - h_3_alpha;
	return solve(X_tau, flux);

}

fc2d_hps_vector<double> merge_w_alpha_refined(fc2d_hps_matrix<double>& X_tau, fc2d_hps_vector<double>& h_3_alpha, fc2d_hps_vector<double>& h_3_beta, fc2d_hps_matrix<double>& L21) {
    fc2d_hps_vector<double> temp = L21 * h_3_alpha;
    temp = h_3_beta - temp;
    return solve(X_tau, temp);
}

fc2d_hps_vector<double> merge_w_beta_refined(fc2d_hps_matrix<double>& X_tau, fc2d_hps_vector<double>& h_3_alpha, fc2d_hps_vector<double>& h_3_beta, fc2d_hps_matrix<double>& L21) {
    fc2d_hps_vector<double> temp = L21 * h_3_beta;
    temp = temp - h_3_alpha;
    return solve(X_tau, temp);
}

fc2d_hps_vector<double> merge_h(fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_vector<double>& w_tau, fc2d_hps_vector<double>& h_1_alpha, fc2d_hps_vector<double>& h_2_beta) {
	fc2d_hps_matrix<double> H(T_13_alpha.rows + T_23_beta.rows, T_13_alpha.cols);
	H.intract(0, 0, T_13_alpha);
	H.intract(T_13_alpha.rows, 0, T_23_beta);
	fc2d_hps_vector<double> h = H * w_tau;

	fc2d_hps_vector<double> temp(h_1_alpha.size() + h_2_beta.size());
	temp.intract(0, h_1_alpha);
	temp.intract(h_1_alpha.size(), h_2_beta);
	h = h + temp;

	return h;
}

fc2d_hps_vector<double> merge_h_alpha_refined(fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_vector<double>& w_tau, fc2d_hps_vector<double>& h_1_alpha, fc2d_hps_vector<double>& h_2_beta, fc2d_hps_matrix<double>& L12) {
    fc2d_hps_matrix<double> temp = T_13_alpha * L12;
    fc2d_hps_matrix<double> H(temp.rows + T_23_beta.rows, temp.cols);
    H.intract(0, 0, temp);
    H.intract(temp.rows, 0, T_23_beta);
    fc2d_hps_vector<double> h = H * w_tau;

    fc2d_hps_vector<double> temp2(h_1_alpha.size() + h_2_beta.size());
	temp2.intract(0, h_1_alpha);
	temp2.intract(h_1_alpha.size(), h_2_beta);
	h = h + temp2;

    return h;
}

fc2d_hps_vector<double> merge_h_beta_refined(fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_vector<double>& w_tau, fc2d_hps_vector<double>& h_1_alpha, fc2d_hps_vector<double>& h_2_beta, fc2d_hps_matrix<double>& L12) {
    fc2d_hps_matrix<double> temp = T_23_beta * L12;
    fc2d_hps_matrix<double> H(T_13_alpha.rows + temp.rows, T_13_alpha.cols);
    H.intract(0, 0, T_13_alpha);
    H.intract(T_13_alpha.rows, 0, temp);
    fc2d_hps_vector<double> h = H * w_tau;

    fc2d_hps_vector<double> temp2(h_1_alpha.size() + h_2_beta.size());
	temp2.intract(0, h_1_alpha);
	temp2.intract(h_1_alpha.size(), h_2_beta);
	h = h + temp2;

    return h;
}

fc2d_hps_patch merge_horizontal_upwards(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

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
    fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I3_alpha);
    fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I3_alpha);

    fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(index_sets.I2, index_sets.I3_beta);
    fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I3_beta);

    fc2d_hps_vector<double> h_1_alpha = alpha.h.from_index_set(index_sets.I1);
	fc2d_hps_vector<double> h_2_beta = beta.h.from_index_set(index_sets.I2);
	fc2d_hps_vector<double> h_3_alpha = alpha.h.from_index_set(index_sets.I3_alpha);
	fc2d_hps_vector<double> h_3_beta = beta.h.from_index_set(index_sets.I3_beta);

    // Perform linear algebra
    fc2d_hps_matrix<double> X_tau; // Is this already built from build stage?
	fc2d_hps_vector<double> w_tau;
	fc2d_hps_vector<double> h_tau;
    if (alpha.N_patch_side[WEST] == beta.N_patch_side[EAST]) {
		// Uniform alpha and beta
        X_tau = merge_X(T_33_alpha, T_33_beta);
		w_tau = merge_w(X_tau, h_3_alpha, h_3_beta);
		h_tau = merge_h(T_13_alpha, T_23_beta, w_tau, h_1_alpha, h_2_beta);
	}
	else if (2 * alpha.N_patch_side[WEST] == beta.N_patch_side[EAST]) {
		// Alpha = fine, beta = coarse
		throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else if (alpha.N_patch_side[WEST] == 2 * beta.N_patch_side[EAST]) {
		// Alpha = coarse, beta = fine
		throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else {
		throw std::invalid_argument("[fc2d_hps_merge::merge_horizontal] Size mismatch between alpha and beta.");
	}

    // Reorder h
    std::vector<int> pi_H = {0, 3, 1, 4, 2, 5};
    std::vector<int> C_S_tau = {N_W_alpha, N_S_alpha, N_N_alpha, N_E_beta, N_S_beta, N_N_beta};
    h_tau = h_tau.block_permute(pi_H, C_S_tau);

    // Create patch
    fc2d_hps_patch merged;
    merged.N_patch_side = {
		alpha.N_patch_side[WEST],
		beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH] + beta.N_patch_side[SOUTH],
		alpha.N_patch_side[NORTH] + beta.N_patch_side[NORTH]
	};
	merged.N_cells_leaf = alpha.N_cells_leaf;
    merged.w = w_tau;
	merged.h = h_tau;

    // Set alpha to hold horizontal merge data
    alpha.w_prime = w_tau; 

    return merged;

}

fc2d_hps_patch merge_vertical_upwards(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

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
    fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(index_sets.I1, index_sets.I3_alpha);
    fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(index_sets.I3_alpha, index_sets.I3_alpha);

    fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(index_sets.I2, index_sets.I3_beta);
    fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(index_sets.I3_beta, index_sets.I3_beta);

	fc2d_hps_vector<double> h_1_alpha = alpha.h.from_index_set(index_sets.I1);
	fc2d_hps_vector<double> h_2_beta = beta.h.from_index_set(index_sets.I2);
	fc2d_hps_vector<double> h_3_alpha = alpha.h.from_index_set(index_sets.I3_alpha);
	fc2d_hps_vector<double> h_3_beta = beta.h.from_index_set(index_sets.I3_beta);

    // Perform linear algebra
    fc2d_hps_matrix<double> X_tau; // Is this already built from build stage?
    fc2d_hps_vector<double> w_tau;
	fc2d_hps_vector<double> h_tau;
	// Begin cases
	if (alpha.N_patch_side[NORTH] == beta.N_patch_side[SOUTH]) {
		// Uniform alpha and beta
        X_tau = merge_X(T_33_alpha, T_33_beta);
		w_tau = merge_w(X_tau, h_3_alpha, h_3_beta);
		h_tau = merge_h(T_13_alpha, T_23_beta, w_tau, h_1_alpha, h_2_beta);
	}
	else if (2 * alpha.N_patch_side[NORTH] == beta.N_patch_side[SOUTH]) {
		// Alpha = fine, beta = coarse
		throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else if (alpha.N_patch_side[NORTH] == 2 * beta.N_patch_side[SOUTH]) {
		// Alpha = coarse, beta = fine
		throw std::logic_error("[fc2d_hps_merge::merge_horizontal] Not implemented!");
	}
	else {
		throw std::invalid_argument("[fc2d_hps_merge::merge_horizontal] Size mismatch between alpha and beta.");
	}

    // Reorder h
    std::vector<int> pi_V = {0, 3, 1, 4, 2, 5};
    std::vector<int> C_S_tau = {N_W_alpha, N_E_alpha, N_S_alpha, N_W_beta, N_E_beta, N_N_beta};
    h_tau = h_tau.block_permute(pi_V, C_S_tau);

    // Create patch
    fc2d_hps_patch merged;
    merged.N_patch_side = {
		alpha.N_patch_side[WEST] + beta.N_patch_side[WEST],
		alpha.N_patch_side[EAST] + beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH],
		beta.N_patch_side[NORTH]
	};
	merged.N_cells_leaf = alpha.N_cells_leaf;
    merged.w = w_tau;
	merged.h = h_tau;

    return merged;

}

void visit_set_particular_data_leaves(fc2d_hps_patch& patch) {

    if (patch.is_leaf) {

        // Create patch solver
        fc2d_hps_FISHPACK_solver FISHPACK_solver;

        fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
        fc2d_hps_options* hps_opt = fc2d_hps_get_options(glob);

        // Get RHS data and set to f
        fclaw2d_domain_t* domain = glob->domain;
        fclaw2d_patch_t* fc_patch = &(domain->blocks->patches[patch.ID]);
        int mfields, meqn;
        double* rhs;
        double* q;
        fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
        fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
        patch.f = fc2d_hps_vector<double>(patch.grid.Nx * patch.grid.Ny);
        for (int i = 0; i < patch.grid.Nx; i++) {
            for (int j = 0; j < patch.grid.Ny; j++) {
                int idx = j + i*patch.grid.Ny;
                int idx_T = i + j*patch.grid.Nx;
                patch.f[idx] = rhs[idx_T];
            }
        }

        // Set Neumann data for particular solution
        fc2d_hps_vector<double> g_zero(2*patch.grid.Nx + 2*patch.grid.Ny, 0);
        patch.h = FISHPACK_solver.dtn(patch.grid, g_zero, patch.f);

        // patch.print_info();

    }
}

void coarsen_patch_upwards(fc2d_hps_patch& fine_patch) {

	// Build L21
	int N_fine = fine_patch.N_cells_leaf * fine_patch.N_patch_side[WEST];
	int N_coarse = N_fine / 2;
	fc2d_hps_matrix<double> L21_side = build_L21(N_coarse, N_fine);
	std::vector<fc2d_hps_matrix<double>> L21_diagonals = {L21_side, L21_side, L21_side, L21_side};
	fc2d_hps_matrix<double> L21_patch = block_diag(L21_diagonals);

	// Particular Neumann data
	// printf("HERE1\n");
	fine_patch.coarsened->h = L21_patch * fine_patch.h;

	// Particular solution data
	// printf("HERE2\n");
	fine_patch.coarsened->w = L21_side * fine_patch.w;

}

void merge_4to1_upwards(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

    // Forward declarations
	fc2d_hps_patch alpha_prime;
	fc2d_hps_patch beta_prime;
	fc2d_hps_patch tau;

	// Check for adaptivity
	fc2d_hps_patch* alpha = &child0;
	fc2d_hps_patch* beta = &child1;
	fc2d_hps_patch* gamma = &child2;
	fc2d_hps_patch* omega = &child3;
	std::vector<int> tags = tag_patch_coarsen(parent, child0, child1, child2, child3);
	// printf("TAGS UPWARDS: ");
	// for (auto& t : tags) printf("%i ", t);
	// printf("\n");
	// if (tags[0]) coarsen_patch_upwards(child0);
	// if (tags[1]) coarsen_patch_upwards(child1);
	// if (tags[2]) coarsen_patch_upwards(child2);
	// if (tags[3]) coarsen_patch_upwards(child3);

	while (tags[0]-- > 0) {
		coarsen_patch_upwards(*alpha);
		alpha = alpha->coarsened;
	}
	while (tags[1]-- > 0) {
		coarsen_patch_upwards(*beta);
		beta = beta->coarsened;
	}
	while (tags[2]-- > 0) {
		coarsen_patch_upwards(*gamma);
		gamma = gamma->coarsened;
	}
	while (tags[3]-- > 0) {
		coarsen_patch_upwards(*omega);
		omega = omega->coarsened;
	}

	// Horizontal merge
	// printf("HERE3\n");
	alpha_prime = merge_horizontal_upwards(*alpha, *beta);
    alpha_prime.T = alpha->T_prime;

	// printf("HERE4\n");
	beta_prime = merge_horizontal_upwards(*gamma, *omega);
    beta_prime.T = gamma->T_prime;

	// Vertical merge
	// printf("HERE5\n");
	tau = merge_vertical_upwards(alpha_prime, beta_prime);
	// printf("HERE6\n");
	
	// Copy only necessary data to parent
    parent.w = tau.w;
    parent.h = tau.h;

}

void fc2d_hps_upwards(fclaw2d_global_t* glob) {
    fclaw_global_essentialf("Begin HPS upwards pass\n");

    // Get quadtree
    fc2d_hps_quadtree<fc2d_hps_patch>* quadtree = fc2d_hps_quadtree<fc2d_hps_patch>::get_instance();

    // Traverse inorder to set non-homogeneous RHS and BC data into solution vector
    quadtree->traverse_postorder(visit_set_particular_data_leaves);

    // Merge upwards to merge non-homogeneous data
    quadtree->merge(merge_4to1_upwards);

    fclaw_global_essentialf("End HPS upwards pass\n");
}