#ifndef FC2D_HPS_POISSONPROBLEM_HPP
#define FC2D_HPS_POISSONPROBLEM_HPP

#include <iostream>
#include <cmath>

enum PROBLEM_TYPE {
	CONSTANT,
	LINEAR,
	LAPLACE1,
	QUAD,
	POLY,
	TRIG,
	GAUSSIAN
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

class PoissonProblem {

public:

	double x_lower, x_upper, y_lower, y_upper;
	PoissonProblem(double x_lower, double x_upper, double y_lower, double y_upper);
	virtual double u(double x, double y) = 0;
	virtual double f(double x, double y) = 0;
	virtual double dudx(double x, double y) = 0;
	virtual double dudy(double x, double y) = 0;

private:

	PoissonProblem(); // Disable default constructor so user has to set problem domain

};

class ConstantLaplaceProblem : PoissonProblem {

public:

	double C;
	ConstantLaplaceProblem(double C, double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double f(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

#endif // FC2D_HPS_POISSONPROBLEM_HPP