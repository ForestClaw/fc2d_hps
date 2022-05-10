#ifndef FC2D_HPP_PATCHSOLVER_HPP
#define FC2D_HPP_PATCHSOLVER_HPP

#include <iostream>
#include <string>
#include <cmath>
#include "fc2d_hps_vector.hpp"
#include "fc2d_hps_matrix.hpp"
#include "fc2d_hps_patchgrid.hpp"

class fc2d_hps_patchsolver {

public:

	virtual std::string get_moniker() = 0;

	virtual fc2d_hps_vector<double> solve(
		fc2d_hps_patchgrid grid,
		fc2d_hps_vector<double> dirichlet_data,
		fc2d_hps_vector<double> rhs_data
	) = 0;

	virtual fc2d_hps_vector<double> dtn(
		fc2d_hps_patchgrid grid,
		fc2d_hps_vector<double> dirichlet_data,
		fc2d_hps_vector<double> rhs_data
	) = 0;

	virtual fc2d_hps_matrix<double> build_dtn(
		fc2d_hps_patchgrid grid
	) = 0;


};

/*
 *	FISHPACK90 Patch Solver
 */

extern "C" {
	void hstcrt_(double* A, double* B, int* M, int* MBDCND, double* BDA, double* BDB, double* C, double* D, int* N, int* NBDCND, double* BDC, double* BDD, double* ELMBDA, double* F, int* IDIMF, double* PERTRB, int* IERROR);
	void hstcrtt_(double* A, double* B, int* M, int* MBDCND, double* BDA, double* BDB, double* C, double* D, int* N, int* NBDCND, double* BDC, double* BDD, double* ELMBDA, double* F, int* IDIMF, double* PERTRB, int* IERROR, double* W);
}

class fc2d_hps_FISHPACK_solver : public fc2d_hps_patchsolver {

public:

	fc2d_hps_FISHPACK_solver();

	std::string get_moniker();
	
	fc2d_hps_vector<double> solve(
		fc2d_hps_patchgrid grid,
		fc2d_hps_vector<double> dirichlet_data,
		fc2d_hps_vector<double> rhs_data
	);

	fc2d_hps_vector<double> dtn(
		fc2d_hps_patchgrid grid,
		fc2d_hps_vector<double> dirichlet_data,
		fc2d_hps_vector<double> rhs_data
	);

	fc2d_hps_matrix<double> build_dtn(
		fc2d_hps_patchgrid grid
	);

};

#endif // FC2D_HPP_PATCHSOLVER_HPP