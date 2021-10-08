#ifndef FC2D_HPS_POISSONPROBLEM_HPP
#define FC2D_HPS_POISSONPROBLEM_HPP

#include <iostream>
#include <cmath>

enum PROBLEM_TYPE {
	CONSTANT,
	LINEAR,
	LAPLACE,
	QUAD,
	POLY,
	TRIG
};

enum BOUNDARY_SIDE {
	WEST,
	EAST,
	SOUTH,
	NORTH
};

class fc2d_hps_poisson_problem {

public:

	int ID;
	double x_lower;
	double x_upper;
	double y_lower;
	double y_upper;

	fc2d_hps_poisson_problem(int problem_ID, double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double f(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

#endif // FC2D_HPS_POISSONPROBLEM_HPP