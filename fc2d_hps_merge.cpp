#include <fc2d_hps_merge.hpp>

fc2d_hps_matrix<double> merge_operation_S(fc2d_hps_matrix<double> T_alpha_33, fc2d_hps_matrix<double> T_beta_33, fc2d_hps_matrix<double> T_alpha_31, fc2d_hps_matrix<double> T_beta_32) {

	fc2d_hps_matrix<double> S_RHS(T_alpha_31.rows, T_alpha_31.cols + T_beta_32.cols);
	// T_alpha_31.negate();
	S_RHS.intract(0, 0, T_alpha_31);
	S_RHS.intract(0, T_alpha_31.cols, T_beta_32);
	fc2d_hps_matrix<double> S = solve(T_alpha_33 - T_beta_33, S_RHS);
	return S;

}

fc2d_hps_matrix<double> merge_operation_T(fc2d_hps_matrix<double> T_alpha_11, fc2d_hps_matrix<double> T_beta_22, fc2d_hps_matrix<double> T_alpha_13, fc2d_hps_matrix<double> T_beta_23, fc2d_hps_matrix<double> S) {

	fc2d_hps_matrix<double> T(T_alpha_11.rows + T_beta_22.rows, T_alpha_11.cols + T_beta_22.cols, 0);
	T.intract(0, 0, T_alpha_11);
	T.intract(T_alpha_11.rows, T_alpha_11.cols, T_beta_22);
	fc2d_hps_matrix<double> T_RHS(T_alpha_13.rows + T_beta_23.rows, T_alpha_13.cols);
	T_RHS.intract(0, 0, T_alpha_13);
	T_RHS.intract(T_alpha_13.rows, 0, T_beta_23);
	T = T + T_RHS*S;
	return T;

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
	// fc2d_hps_patch lower_merged();
	// fc2d_hps_patch upper_merged();

	fc2d_hps_matrix<double> T_alpha = child0.T; // TODO: Write copy constructor for fc2d_hps_matrix<T>
	fc2d_hps_matrix<double> T_beta = child1.T;

	// ...

	// Vertical merge
	
	

}