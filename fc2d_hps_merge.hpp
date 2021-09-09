#ifndef FC2D_HPS_MERGE_HPP
#define FC2D_HPS_MERGE_HPP

#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"

fc2d_hps_matrix<double> merge_operation_S(
	fc2d_hps_matrix<double> T_alpha_33,
	fc2d_hps_matrix<double> T_beta_33,
	fc2d_hps_matrix<double> T_alpha_31,
	fc2d_hps_matrix<double> T_beta_32
);

fc2d_hps_matrix<double> merge_operation_T(
	fc2d_hps_matrix<double> T_alpha_11,
	fc2d_hps_matrix<double> T_beta_22,
	fc2d_hps_matrix<double> T_alpha_13,
	fc2d_hps_matrix<double> T_beta_23,
	fc2d_hps_matrix<double> S
);

#endif // FC2D_HPS_MERGE_HPP