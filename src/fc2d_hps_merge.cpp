#include "fc2d_hps_merge.hpp"

fc2d_hps_patch merge_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

	// Build index vectors
	int N_points_leaf_side = alpha.grid.Nx; // TODO: Need to make sure this is in fact for a leaf patch
	std::vector<int> I_W_alpha(alpha.N_patch_side[WEST] * N_points_leaf_side);
	std::vector<int> I_E_alpha(alpha.N_patch_side[EAST] * N_points_leaf_side);
	std::vector<int> I_S_alpha(alpha.N_patch_side[SOUTH] * N_points_leaf_side);
	std::vector<int> I_N_alpha(alpha.N_patch_side[NORTH] * N_points_leaf_side);
	std::vector<int> I_W_beta(beta.N_patch_side[WEST] * N_points_leaf_side);
	std::vector<int> I_E_beta(beta.N_patch_side[EAST] * N_points_leaf_side);
	std::vector<int> I_S_beta(beta.N_patch_side[SOUTH] * N_points_leaf_side);
	std::vector<int> I_N_beta(beta.N_patch_side[NORTH] * N_points_leaf_side);
	std::vector<int> I_1(0);
	std::vector<int> I_2(0);
	std::vector<int> I_3_alpha(0);
	std::vector<int> I_3_beta(0);

	std::iota(I_W_alpha.begin(), I_W_alpha.end(), 0);
	std::iota(I_E_alpha.begin(), I_E_alpha.end(), alpha.N_patch_side[WEST]);
	std::iota(I_S_alpha.begin(), I_S_alpha.end(), alpha.N_patch_side[WEST] + alpha.N_patch_side[EAST]);
	std::iota(I_N_alpha.begin(), I_N_alpha.end(), alpha.N_patch_side[WEST] + alpha.N_patch_side[EAST] + alpha.N_patch_side[SOUTH]);
	std::iota(I_W_beta.begin(), I_W_beta.end(), 0);
	std::iota(I_E_beta.begin(), I_E_beta.end(), beta.N_patch_side[WEST]);
	std::iota(I_S_beta.begin(), I_S_beta.end(), beta.N_patch_side[WEST] + beta.N_patch_side[EAST]);
	std::iota(I_N_beta.begin(), I_N_beta.end(), beta.N_patch_side[WEST] + beta.N_patch_side[EAST] + alpha.N_patch_side[SOUTH]);

	I_1.insert(I_1.end(), I_W_alpha.begin(), I_W_alpha.end());
	I_1.insert(I_1.end(), I_S_alpha.begin(), I_S_alpha.end());
	I_1.insert(I_1.end(), I_N_alpha.begin(), I_N_alpha.end());
	I_2.insert(I_2.end(), I_E_beta.begin(), I_E_beta.end());
	I_2.insert(I_2.end(), I_S_beta.begin(), I_S_beta.end());
	I_2.insert(I_2.end(), I_N_beta.begin(), I_N_beta.end());
	I_3_alpha.insert(I_3_alpha.end(), I_E_alpha.begin(), I_E_alpha.end());
	I_3_beta.insert(I_3_beta.end(), I_W_beta.begin(), I_W_beta.end());

	// Extract blocks
	fc2d_hps_matrix<double> T_11_alpha = alpha.T.from_index_set(I_1, I_1);
	fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(I_1, I_3_alpha);
	fc2d_hps_matrix<double> T_31_alpha = alpha.T.from_index_set(I_3_alpha, I_1);
	fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(I_3_alpha, I_3_alpha);

	fc2d_hps_matrix<double> T_22_beta = beta.T.from_index_set(I_2, I_2);
	fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(I_2, I_3_beta);
	fc2d_hps_matrix<double> T_32_beta = beta.T.from_index_set(I_3_beta, I_2);
	fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(I_3_beta, I_3_beta);

	// Perform merge linear algebra
	// @TODO

	// Create new merged patch
	fc2d_hps_patch merged;
	merged.level = alpha.level;
	merged.is_leaf = false;
	merged.N_patch_side = {
		alpha.N_patch_side[WEST],
		beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH] + beta.N_patch_side[SOUTH],
		alpha.N_patch_side[NORTH] + beta.N_patch_side[NORTH]
	};
	// merged.X = X_tau;
	// merged.H = H_tau;
	// merged.S = S_tau;
	// merged.T = T_tau;
	return merged;
	
}

fc2d_hps_patch merge_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

	// Build index vectors
	int N_points_leaf_side = alpha.grid.Nx; // TODO: Need to make sure this is in fact for a leaf patch
	std::vector<int> I_W_alpha(alpha.N_patch_side[WEST] * N_points_leaf_side);
	std::vector<int> I_E_alpha(alpha.N_patch_side[EAST] * N_points_leaf_side);
	std::vector<int> I_S_alpha(alpha.N_patch_side[SOUTH] * N_points_leaf_side);
	std::vector<int> I_N_alpha(alpha.N_patch_side[NORTH] * N_points_leaf_side);
	std::vector<int> I_W_beta(beta.N_patch_side[WEST] * N_points_leaf_side);
	std::vector<int> I_E_beta(beta.N_patch_side[EAST] * N_points_leaf_side);
	std::vector<int> I_S_beta(beta.N_patch_side[SOUTH] * N_points_leaf_side);
	std::vector<int> I_N_beta(beta.N_patch_side[NORTH] * N_points_leaf_side);
	std::vector<int> I_1(0);
	std::vector<int> I_2(0);
	std::vector<int> I_3_alpha(0);
	std::vector<int> I_3_beta(0);

	std::iota(I_W_alpha.begin(), I_W_alpha.end(), 0);
	std::iota(I_E_alpha.begin(), I_E_alpha.end(), alpha.N_patch_side[WEST]);
	std::iota(I_S_alpha.begin(), I_S_alpha.end(), alpha.N_patch_side[WEST] + alpha.N_patch_side[EAST]);
	std::iota(I_N_alpha.begin(), I_N_alpha.end(), alpha.N_patch_side[WEST] + alpha.N_patch_side[EAST] + alpha.N_patch_side[SOUTH]);
	std::iota(I_W_beta.begin(), I_W_beta.end(), 0);
	std::iota(I_E_beta.begin(), I_E_beta.end(), beta.N_patch_side[WEST]);
	std::iota(I_S_beta.begin(), I_S_beta.end(), beta.N_patch_side[WEST] + beta.N_patch_side[EAST]);
	std::iota(I_N_beta.begin(), I_N_beta.end(), beta.N_patch_side[WEST] + beta.N_patch_side[EAST] + alpha.N_patch_side[SOUTH]);

	I_1.insert(I_1.end(), I_W_alpha.begin(), I_W_alpha.end());
	I_1.insert(I_1.end(), I_E_alpha.begin(), I_E_alpha.end());
	I_1.insert(I_1.end(), I_S_alpha.begin(), I_S_alpha.end());
	I_2.insert(I_2.end(), I_W_beta.begin(), I_W_beta.end());
	I_2.insert(I_2.end(), I_E_beta.begin(), I_E_beta.end());
	I_2.insert(I_2.end(), I_N_beta.begin(), I_N_beta.end());
	I_3_alpha.insert(I_3_alpha.end(), I_N_alpha.begin(), I_N_alpha.end());
	I_3_beta.insert(I_3_beta.end(), I_S_beta.begin(), I_S_beta.end());

	// Extract blocks
	fc2d_hps_matrix<double> T_11_alpha = alpha.T.from_index_set(I_1, I_1);
	fc2d_hps_matrix<double> T_13_alpha = alpha.T.from_index_set(I_1, I_3_alpha);
	fc2d_hps_matrix<double> T_31_alpha = alpha.T.from_index_set(I_3_alpha, I_1);
	fc2d_hps_matrix<double> T_33_alpha = alpha.T.from_index_set(I_3_alpha, I_3_alpha);

	fc2d_hps_matrix<double> T_22_beta = beta.T.from_index_set(I_2, I_2);
	fc2d_hps_matrix<double> T_23_beta = beta.T.from_index_set(I_2, I_3_beta);
	fc2d_hps_matrix<double> T_32_beta = beta.T.from_index_set(I_3_beta, I_2);
	fc2d_hps_matrix<double> T_33_beta = beta.T.from_index_set(I_3_beta, I_3_beta);

	// Perform merge linear algebra
	// @TODO

	// Create new merged patch
	fc2d_hps_patch merged;
	merged.level = alpha.level - 1;
	merged.is_leaf = false;
	merged.N_patch_side = {
		alpha.N_patch_side[WEST] + beta.N_patch_side[WEST],
		alpha.N_patch_side[EAST] + beta.N_patch_side[EAST],
		alpha.N_patch_side[SOUTH],
		beta.N_patch_side[NORTH]
	};
	// merged.X = X_tau;
	// merged.H = H_tau;
	// merged.S = S_tau;
	// merged.T = T_tau;
	return merged;

}

void merge_4to1(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

	// Error checks
	//    Check if all children are leaves
	if (child0.is_leaf == false) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child0` is not a leaf"); }
	if (child1.is_leaf == false) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child1` is not a leaf"); }
	if (child2.is_leaf == false) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child2` is not a leaf"); }
	if (child3.is_leaf == false) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child3` is not a leaf"); }

	//    Check if children data are set
	if (child0.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child0.T.size()` is 0; it shouldn't be..."); }
	if (child1.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child1.T.size()` is 0; it shouldn't be..."); }
	if (child2.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child2.T.size()` is 0; it shouldn't be..."); }
	if (child3.T.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child3.T.size()` is 0; it shouldn't be..."); }
	if (child0.S.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child0.S.size()` is 0; it shouldn't be..."); }
	if (child1.S.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child1.S.size()` is 0; it shouldn't be..."); }
	if (child2.S.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child2.S.size()` is 0; it shouldn't be..."); }
	if (child3.S.size() == 0) { throw std::invalid_argument("[fc2d_hps_merge merge_4to1] `child3.S.size()` is 0; it shouldn't be..."); }

	// Horizontal merge
	fc2d_hps_patch alpha_prime = merge_horizontal(child0, child1);
	fc2d_hps_patch beta_prime = merge_horizontal(child2, child3);

	// Vertical merge
	parent = merge_vertical(alpha_prime, beta_prime);

}