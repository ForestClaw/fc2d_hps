#include <Structures/fc2d_hps_poissonproblem.hpp>

// ==================== Base Poisson Class ====================
fc2d_hps_poisson_problem::fc2d_hps_poisson_problem(double x_lower, double x_upper, double y_lower, double y_upper) :
	x_lower(x_lower), x_upper(x_upper), y_lower(y_lower), y_upper(y_upper)
		{}

// ==================== Base Laplace Class ====================
fc2d_hps_laplace_problem::fc2d_hps_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper) :
	fc2d_hps_poisson_problem(x_lower, x_upper, y_lower, y_upper)
		{}

double fc2d_hps_laplace_problem::f(double x, double y) {
	return 0;
}

// ==================== Constant Laplace Equation ====================
fc2d_hps_constant_laplace_problem::fc2d_hps_constant_laplace_problem(double C, double x_lower, double x_upper, double y_lower, double y_upper) :
	fc2d_hps_laplace_problem(x_lower, x_upper, y_lower, y_upper), C(C)
		{}

double fc2d_hps_constant_laplace_problem::u(double x, double y) {
	return this->C;
}

double fc2d_hps_constant_laplace_problem::dudx(double x, double y) {
	return 0.0;
}

double fc2d_hps_constant_laplace_problem::dudy(double x, double y) {
	return 0.0;
}

// ==================== Linear Laplace Equation ====================
fc2d_hps_linear_laplace_problem::fc2d_hps_linear_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper) :
	fc2d_hps_laplace_problem(x_lower, x_upper, y_lower, y_upper)
		{}

double fc2d_hps_linear_laplace_problem::u(double x, double y) {
	return x + y;
}

double fc2d_hps_linear_laplace_problem::dudx(double x, double y) {
	return 1.0;
}

double fc2d_hps_linear_laplace_problem::dudy(double x, double y) {
	return 1.0;
}

// ==================== Hyperbolic Laplace Equation ====================
fc2d_hps_hyperbolic_laplace_problem::fc2d_hps_hyperbolic_laplace_problem(double x_lower, double x_upper, double y_lower, double y_upper) :
	fc2d_hps_laplace_problem(x_lower, x_upper, y_lower, y_upper)
		{}

double fc2d_hps_hyperbolic_laplace_problem::u(double x, double y) {
	return sin(this->b*x) * sinh(this->b*y);
}

double fc2d_hps_hyperbolic_laplace_problem::dudx(double x, double y) {
	return this->b*cos(this->b*x)*sinh(this->b*y);
}

double fc2d_hps_hyperbolic_laplace_problem::dudy(double x, double y) {
	return this->b*cosh(this->b*y)*sin(this->b*x);
}

// ==================== Quadratic Poisson Equation ====================
fc2d_hps_quadratic_poisson_problem::fc2d_hps_quadratic_poisson_problem(double x_lower, double x_upper, double y_lower, double y_upper) :
	fc2d_hps_poisson_problem(x_lower, x_upper, y_lower, y_upper)
		{}

double fc2d_hps_quadratic_poisson_problem::u(double x, double y) {
	return pow(x, 2) + pow(y, 2) + 2*x*y;
}

double fc2d_hps_quadratic_poisson_problem::f(double x, double y) {
	return 4.0;
}

double fc2d_hps_quadratic_poisson_problem::dudx(double x, double y) {
	return 2.0*x + 2.0*y;
}

double fc2d_hps_quadratic_poisson_problem::dudy(double x, double y) {
	return 2.0*x + 2.0*y;
}