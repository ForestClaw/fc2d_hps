#include <Structures/fc2d_hps_patch.hpp>

fc2d_hps_patch::fc2d_hps_patch() {}

fc2d_hps_patch::fc2d_hps_patch(fc2d_hps_patchgrid patch_grid, int ID, int level, bool is_leaf) :
	ID(ID), level(level), is_leaf(is_leaf), grid(patch_grid) {
		if (is_leaf) {
			N_patch_side[WEST] = 1;
			N_patch_side[EAST] = 1;
			N_patch_side[SOUTH] = 1;
			N_patch_side[NORTH] = 1;
		}
	}

void fc2d_hps_patch::print_info() {
	printf("---------- PATCH INFO ----------\n");
	printf("  ID           = %i\n", this->ID);
	printf("  level        = %i\n", this->level);
	printf("  is_leaf      = %i\n", this->is_leaf);
	printf("  has_coarsened= %i\n", this->has_coarsened);
	printf("  has_finer    = %i\n", this->has_finer);
	printf("  N_cells_leaf = %i\n", this->N_cells_leaf);
	printf("  N_patch_side:\n");
	printf("    WEST  = %i\n", this->N_patch_side[WEST]);
	printf("    EAST  = %i\n", this->N_patch_side[EAST]);
	printf("    SOUTH = %i\n", this->N_patch_side[SOUTH]);
	printf("    NORTH = %i\n", this->N_patch_side[NORTH]);
	printf("  Grid:\n");
	printf("    grid.Nx      = %i\n", this->grid.Nx);
	printf("    grid.Ny      = %i\n", this->grid.Ny);
	printf("    grid.dx      = %11.4e\n", this->grid.dx);
	printf("    grid.dy      = %11.4e\n", this->grid.dy);
	printf("    grid.x_lower = %11.4e\n", this->grid.x_lower);
	printf("    grid.x_upper = %11.4e\n", this->grid.x_upper);
	printf("    grid.y_lower = %11.4e\n", this->grid.y_lower);
	printf("    grid.y_upper = %11.4e\n", this->grid.y_upper);
	printf("  Data Sizes:\n");
	printf("    T       : [%8i, %8i]\n", this->T.rows, this->T.cols);
	printf("    S       : [%8i, %8i]\n", this->S.rows, this->S.cols);
	printf("    S_prime : [%8i, %8i]\n", this->S_prime.rows, this->S_prime.cols);
	printf("    X       : [%8i, %8i]\n", this->X.rows, this->X.cols);
	printf("    u       : [%8i]\n", this->u.size());
	printf("    f       : [%8i]\n", this->f.size());
	printf("    g       : [%8i]\n", this->g.size());
	printf("    h       : [%8i]\n", this->h.size());
	printf("    w       : [%8i]\n", this->w.size());
	printf("    w_prime : [%8i]\n", this->w_prime.size());
	printf("--------------------------------\n");
}

void fc2d_hps_patch::write_internal_data(std::ofstream& file, const fc2d_hps_vector<double>& v, std::string data_name) {
	file << "SCALARS " << data_name << " double" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < this->grid.Nx; i++) {
		for (int j = 0; j < this->grid.Ny; j++) {
			int idx = j + i*this->grid.Ny;
			file << v[idx] << std::endl;
		}
	}
	file << std::endl;
}

void fc2d_hps_patch::to_vtk(std::string filename, std::string comments, std::string outputs) {

	// Check if internal data points are wanted
	if (
		outputs.find('u') != std::string::npos ||
		outputs.find('f') != std::string::npos ||
		outputs.find('w') != std::string::npos
	) {
		// Open file
		std::ofstream file(filename + ".vtk");

		// Write header
		file << "# vtk DataFile Version 2.0" << std::endl;
		file << comments << std::endl;
		file << "ASCII" << std::endl;

		// Write grid
		file << "DATASET STRUCTURED_GRID" << std::endl;
		file << "DIMENSIONS " << this->grid.Nx << " " << this->grid.Ny << " 1" << std::endl;
		file << "POINTS " << this->grid.Nx*this->grid.Ny << " double" << std::endl;
		for (int i = 0; i < this->grid.Nx; i++) {
			for (int j = 0; j < this->grid.Ny; j++) {
				file << this->grid.point(XDIM, i) << " " << this->grid.point(YDIM, j) << " 0" << std::endl;
			}
		}
		file << std::endl;

		// Write point data
		file << "POINT_DATA " << this->grid.Nx*this->grid.Ny << std::endl;

		// Write u data
		if (outputs.find('u') != std::string::npos) {
			write_internal_data(file, this->u, "u_SOLUTION_DATA");
		}

		// Write f data
		if (outputs.find('f') != std::string::npos) {
			write_internal_data(file, this->f, "f_RHS_DATA");
		}

		// Write w data
		if (outputs.find('w') != std::string::npos) {
			write_internal_data(file, this->w, "w_PARTICULAR_DATA");
		}

		// Close file
		file.close();

	}

}

// std::ofstream quadtree_file;
// std::vector<int> number_points(3, 0);

// void visit_setup_vtk(fc2d_hps_patch& patch) {
// 	if (patch.is_leaf) {
// 		number_points[0] += patch.grid.Nx;
// 		number_points[1] += patch.grid.Ny;
// 		number_points[2] = 1;
// 	}
// }

// void visit_write_grid(fc2d_hps_patch& patch) {
// 	if (patch.is_leaf) {
// 		for (int i = 0; i < patch.grid.Nx; i++) {
// 			for (int j = 0; j < patch.grid.Ny; j++) {
// 				quadtree_file << patch.grid.point(XDIM, i) << " " << patch.grid.point(YDIM, j) << " 0" << std::endl;
// 			}
// 		}
// 	}
// }

// void visit_write_data(fc2d_hps_patch& patch) {
// 	if (patch.is_leaf) {
// 		for (int i = 0; i < patch.grid.Nx; i++) {
// 			for (int j = 0; j < patch.grid.Ny; j++) {
// 				int idx = j + i*patch.grid.Ny;
// 				quadtree_file << patch.u[idx] << std::endl;
// 			}
// 		}
// 	}
// }

// void hps_patch_quadtree_to_vtk(fc2d_hps_quadtree<fc2d_hps_patch>& quadtree, std::string filename, std::string comments) {
	
// 	// Create file to write to
// 	quadtree_file = std::ofstream(filename + ".vtk");

// 	// Write header
// 	quadtree_file << "# vtk DataFile Version 2.0" << std::endl;
// 	quadtree_file << comments << std::endl;
// 	quadtree_file << "ASCII" << std::endl;

// 	// Setup values
// 	quadtree.traverse_inorder(visit_setup_vtk);

// 	// Write grid header
// 	int total_points = 1;
// 	for (auto& n : number_points) total_points *= n;
// 	quadtree_file << "DATASET STRUCTURED_GRID" << std::endl;
// 	quadtree_file << "DIMENSIONS " << number_points[0] << " " << number_points[1] << " " << number_points[2] << std::endl;
// 	quadtree_file << "POINTS " << total_points << " double" << std::endl;

// 	// Write grid for each patch
// 	quadtree.traverse_inorder(visit_write_grid);

// 	// Write data header
// 	quadtree_file << "POINT_DATA " << total_points << std::endl;

// 	// Write data for each patch
// 	//    u solution data
// 	quadtree_file << "SCALARS u_SOLUTION_DATA double" << std::endl;
// 	quadtree_file << "LOOKUP_TABLE default" << std::endl;
// 	quadtree.traverse_inorder(visit_write_data);
// 	quadtree_file << std::endl;

// 	// @TODO: Others...

// 	// Close file
// 	quadtree_file.close();
// }