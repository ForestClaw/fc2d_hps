#ifndef FC2D_HPS_MERGE_HPP
#define FC2D_HPS_MERGE_HPP

#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fclaw2d_global.h>
#include <Util/fc2d_hps_options.h>
#include <Structures/fc2d_hps_vector.hpp>
#include <Structures/fc2d_hps_matrix.hpp>
#include <Structures/fc2d_hps_patch.hpp>
typedef struct index_set {

    std::vector<int> I1;
    std::vector<int> I2;
    std::vector<int> I3_alpha;
    std::vector<int> I3_beta;

} index_set_t;

std::vector<int> fill_range(int start, int end);
fc2d_hps_matrix<double> build_L21(int n_rows, int n_cols);
fc2d_hps_matrix<double> build_L12(int n_rows, int n_cols);
fc2d_hps_matrix<double> merge_X(fc2d_hps_matrix<double>& T_33_alpha, fc2d_hps_matrix<double>& T_33_beta);
fc2d_hps_matrix<double> merge_S(fc2d_hps_matrix<double>& X_tau, fc2d_hps_matrix<double>& T_31_alpha, fc2d_hps_matrix<double>& T_32_beta);
fc2d_hps_matrix<double> merge_T(fc2d_hps_matrix<double>& S_tau, fc2d_hps_matrix<double>& T_11_alpha, fc2d_hps_matrix<double>& T_22_beta, fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta);
fc2d_hps_vector<double> merge_w(fc2d_hps_matrix<double>& X_tau, fc2d_hps_vector<double>& h_3_alpha, fc2d_hps_vector<double>& h_3_beta);
fc2d_hps_vector<double> merge_h(fc2d_hps_matrix<double>& T_13_alpha, fc2d_hps_matrix<double>& T_23_beta, fc2d_hps_vector<double>& w_tau, fc2d_hps_vector<double>& h_1_alpha, fc2d_hps_vector<double>& h_2_beta);
index_set_t make_index_sets_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
fc2d_hps_patch merge_horizontal(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
index_set_t make_index_sets_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
fc2d_hps_patch merge_vertical(fc2d_hps_patch& alpha, fc2d_hps_patch& beta);
std::vector<int> tag_patch_coarsen(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3);
void merge_4to1(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3);

#endif // FC2D_HPS_MERGE_HPP