#ifndef FC2D_HPS_PATCH_HPP
#define FC2D_HPS_PATCH_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"
#include "fc2d_hps_quadtree.hpp"
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

	// Data matrices
	fc2d_hps_matrix<double> T;			// DtN Matrix
	fc2d_hps_matrix<double> S;			// Solution Matrix
	fc2d_hps_matrix<double> S_prime;	// Horizontal Merge Solution Matrix
	fc2d_hps_matrix<double> X;			// Body Load Matrix
	
	// Solution data
	fc2d_hps_vector<double> u;			// Solution Vector
	fc2d_hps_vector<double> g;			// Dirichlet Vector
	fc2d_hps_vector<double> v;			// Neumann Vector

	// Homogeneous data
	fc2d_hps_vector<double> f;			// Poisson Vector

	// Particular data
	fc2d_hps_vector<double> h;			// Particular Neumann Vector
	fc2d_hps_vector<double> w;			// Particular Solution Vector
	fc2d_hps_vector<double> w_prime;	// Horizontal Particular Solution Vector

	fc2d_hps_patch();
	fc2d_hps_patch(fc2d_hps_patchgrid patch_grid, int ID, int level, bool is_leaf);

	// @TODO: Output functions
	void print_info();
	void to_mmio();
	void to_vtk(std::string filename, std::string comments="Output Path Info", std::string outputs="");

private:

	void write_internal_data(std::ofstream& file, const fc2d_hps_vector<double>& v, std::string data_name);

};

extern std::ofstream quadtree_file;
void visit_setup_vtk(fc2d_hps_patch& patch);
void visit_write_grid(fc2d_hps_patch& patch);
void visit_write_data(fc2d_hps_patch& patch);
void hps_patch_quadtree_to_vtk(fc2d_hps_quadtree<fc2d_hps_patch>& quadtree, std::string filename, std::string comments="HPS Patch Quadtree Output");

#endif // FC2D_HPS_PATCH_HPP