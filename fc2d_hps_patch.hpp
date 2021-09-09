#ifndef FC2D_HPS_PATCH_HPP
#define FC2D_HPS_PATCH_HPP

#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"

class fc2d_hps_patch {

public:

	// Metadata
	int ID;
	int level;
	bool is_leaf;

	// TODO: Data structure or pointer to store patch grid information

	// HPS Data
	fc2d_hps_matrix<double> T;	// DtN Matrix
	fc2d_hps_matrix<double> S;	// Solution Matrix
	fc2d_hps_matrix<double> F;	// Body Load Matrix (placeholder)
	fc2d_hps_matrix<double> H;	// Flux Matrix (placeholder)
	fc2d_hps_vector<double> u;	// Solution Vector
	fc2d_hps_vector<double> f;	// Poisson Vector
	fc2d_hps_vector<double> g;	// Dirichlet Vector
	fc2d_hps_vector<double> h;	// Neumann Vector
	fc2d_hps_vector<double> w;	// Particular Solution Vector

	fc2d_hps_patch();

};

#endif // FC2D_HPS_PATCH_HPP