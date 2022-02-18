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

/**
 * 
 * 
 */
class fc2d_hps_poisson_problem {

public:

	double x_lower, x_upper, y_lower, y_upper;
	fc2d_hps_poisson_problem(double x_lower, double x_upper, double y_lower, double y_upper);
	virtual double u(double x, double y) = 0;
	virtual double f(double x, double y) = 0;
	virtual double dudx(double x, double y) = 0;
	virtual double dudy(double x, double y) = 0;

private:

	fc2d_hps_poisson_problem(); // Disable default constructor so user has to set problem domain

};

/**
 * 
 * 
 */
class fc2d_hps_laplace_problem : public fc2d_hps_poisson_problem {

public:

	fc2d_hps_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper);
	virtual double u(double x, double y) = 0;
	virtual double dudx(double x, double y) = 0;
	virtual double dudy(double x, double y) = 0;
	double f(double x, double y);

};

/**
 * 
 * 
 */
class fc2d_hps_constant_laplace_problem : public fc2d_hps_laplace_problem {

public:

	double C;
	fc2d_hps_constant_laplace_problem(double C, double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

/**
 * 
 * 
 */
class fc2d_hps_linear_laplace_problem : public fc2d_hps_laplace_problem {

public:

	fc2d_hps_linear_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

/**
 * 
 * 
 */
class fc2d_hps_hyperbolic_laplace_problem : public fc2d_hps_laplace_problem {

private:

	double b = (2.0 / 3.0) * M_PI;

public:

	fc2d_hps_hyperbolic_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

class fc2d_hps_quadratic_poisson_problem : public fc2d_hps_poisson_problem {

public:

	fc2d_hps_quadratic_poisson_problem(double x_lower, double x_upper, double y_lower, double y_upper);
	double u(double x, double y);
	double f(double x, double y);
	double dudx(double x, double y);
	double dudy(double x, double y);

};

#endif // FC2D_HPS_POISSONPROBLEM_HPP