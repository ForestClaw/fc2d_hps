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