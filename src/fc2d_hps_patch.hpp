#ifndef FC2D_HPS_PATCH_HPP
#define FC2D_HPS_PATCH_HPP

#include <vector>
#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"
#include "fc2d_hps_patchgrid.hpp"
#include "fc2d_hps_poissonproblem.hpp" // For sides enum

// enum PATCH_SIDES {
// 	WEST,
// 	EAST,
// 	SOUTH,
// 	NORTH
// };

class fc2d_hps_patch {

public:

	// Metadata
	int ID;											// Patch's global ID
	int level;										// Level in tree
	bool is_leaf;									// Flag for if patch is a leaf
	int N_cells_leaf;								// Storage for number of cells on leaf patch side
	std::vector<int> N_patch_side = {0, 0, 0, 0};	// To keep track of patch's side based on children

	// Patch grid information
	fc2d_hps_patchgrid grid;	// Grid information

	// User data
	void* user;

	// HPS Data
	fc2d_hps_matrix<double> T;			// DtN Matrix
	fc2d_hps_matrix<double> S;			// Solution Matrix
	fc2d_hps_matrix<double> S_prime;	// Horizontal Merge Solution Matrix
	fc2d_hps_matrix<double> X;			// Body Load Matrix
	fc2d_hps_matrix<double> H;			// Flux Matrix (placeholder)
	fc2d_hps_vector<double> u;			// Solution Vector
	fc2d_hps_vector<double> f;			// Poisson Vector
	fc2d_hps_vector<double> g;			// Dirichlet Vector
	fc2d_hps_vector<double> h;			// Neumann Vector
	fc2d_hps_vector<double> w;			// Particular Solution Vector
	fc2d_hps_vector<double> w_prime;	// Horizontal Particular Solution Vector

	fc2d_hps_patch();
	fc2d_hps_patch(fc2d_hps_patchgrid patch_grid, int ID, int level, bool is_leaf);

	// @TODO: Output functions
	void print_info();

};

#endif // FC2D_HPS_PATCH_HPP