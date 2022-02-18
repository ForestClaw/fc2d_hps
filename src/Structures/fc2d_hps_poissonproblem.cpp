#include <Structures/fc2d_hps_poissonproblem.hpp>
/*
// fc2d_hps_poisson_problem::fc2d_hps_poisson_problem(int problem_ID, double x_lower, double x_upper, double y_lower, double y_upper) :
// 	ID(problem_ID), x_lower(x_lower), x_upper(x_upper), y_lower(y_lower), y_upper(y_upper)
// 		{}

// double fc2d_hps_poisson_problem::u(double x, double y) {
// 	double b = (2.0/3.0) * M_PI;
// 	double A = 1.0;
// 	double x0 = (this->x_lower + this->x_upper) / 2.0;
// 	double y0 = (this->y_lower + this->y_upper) / 2.0;
// 	double sigma_x = 0.3;
// 	double sigma_y = 0.3;
// 	double x_part = pow(x - x0, 2) / (2.0*pow(sigma_x, 2));
// 	double y_part = pow(y - y0, 2) / (2.0*pow(sigma_y, 2));
// 	switch(ID) {
// 		case CONSTANT:
// 			return 1.0;
// 		case LINEAR:
// 			return x;
// 		case LAPLACE1:
// 			return sin(b*x) * sinh(b*y);
// 		case QUAD:
// 			return x*x + y*y + 2*x*y;
// 		case POLY:
// 			return y*y*x*x*x*x;
// 		case TRIG:
// 			return sin(2*M_PI*x) * sin(2*M_PI*y);
// 		case GAUSSIAN:
// 			return A*exp(-(x_part + y_part));
// 		default:
// 			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
// 	}
// }

// double fc2d_hps_poisson_problem::f(double x, double y) {
// 	double A = 1.0;
// 	double x0 = (this->x_lower + this->x_upper) / 2.0;
// 	double y0 = (this->y_lower + this->y_upper) / 2.0;
// 	double sigma_x = 0.3;
// 	double sigma_y = 0.3;
// 	double x_part = pow(x - x0, 2) / (2*pow(sigma_x, 2));
// 	double y_part = pow(y - y0, 2) / (2*pow(sigma_y, 2));
// 	switch(ID) {
// 		case CONSTANT:
// 			return 0.0;
// 		case LINEAR:
// 			return 0.0;
// 		case LAPLACE1:
// 			return 0.0;
// 		case QUAD:
// 			return 4.0;
// 		case POLY:
// 			return 2.0*pow(x, 4) + 12.0*pow(x,2)*pow(y,2);
// 		case TRIG:
// 			return -2.0*pow(2.0*M_PI,2)*u(x,y);
// 		case GAUSSIAN:
// 			return u(x,y) * ((pow(y,2)*pow(sigma_y,4) - 2.0*y*y0*pow(sigma_x,4) + pow(y0,2)*pow(sigma_x,4) - pow(sigma_x,4)*pow(sigma_y,2) + pow(x,2)*pow(sigma_y,4) - 2.0*x*x0*pow(sigma_y,4) + pow(x0,2)*pow(sigma_y,4) - pow(sigma_x,2)*pow(sigma_y,4)) / (pow(sigma_x*sigma_y,4)));
// 		default:
// 			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
// 	}

// }

// double fc2d_hps_poisson_problem::dudx(double x, double y) {
// 	double b = (2.0/3.0)*M_PI;
// 	double x0 = (this->x_lower + this->x_upper) / 2.0;
// 	double y0 = (this->y_lower + this->y_upper) / 2.0;
// 	double sigma_x = 0.3;
// 	double sigma_y = 0.3;
// 	double x_part = pow(x - x0, 2) / (2*pow(sigma_x, 2));
// 	double y_part = pow(y - y0, 2) / (2*pow(sigma_y, 2));
// 	switch(ID) {
// 		case CONSTANT:
// 			return 0.0;
// 		case LINEAR:
// 			return 1.0;
// 		case LAPLACE1:
// 			return b*cos(b*x)*sinh(b*y);
// 		case QUAD:
// 			return 2.0*x + 2.0*y;
// 		case POLY:
// 			return 4.0*y*y*x*x*x;
// 		case TRIG:
// 			return 2.0*M_PI*cos(2.0*M_PI*x) * sin(2.0*M_PI*y);
// 		case GAUSSIAN:
// 			return -u(x,y) * (x - x0) / pow(sigma_x,2);
// 		default:
// 			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
// 	}
// }

// double fc2d_hps_poisson_problem::dudy(double x, double y) {
// 	double b = (2.0/3.0)*M_PI;
// 	double x0 = (this->x_lower + this->x_upper) / 2.0;
// 	double y0 = (this->y_lower + this->y_upper) / 2.0;
// 	double sigma_x = 0.3;
// 	double sigma_y = 0.3;
// 	double x_part = pow(x - x0, 2) / (2*pow(sigma_x, 2));
// 	double y_part = pow(y - y0, 2) / (2*pow(sigma_y, 2));
// 	switch(ID) {
// 		case CONSTANT:
// 			return 0.0;
// 		case LINEAR:
// 			return 0.0;
// 		case LAPLACE1:
// 			return b*cosh(b*y)*sin(b*x);
// 		case QUAD:
// 			return 2.0*y + 2*x;
// 		case POLY:
// 			return 2.0*y*x*x*x*x;
// 		case TRIG:
// 			return sin(2*M_PI*x) * 2.0*M_PI*cos(2*M_PI*y);
// 		case GAUSSIAN:
// 			return -u(x,y) * (y - y0) / pow(sigma_y,2);
// 		default:
// 			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
// 	}
// }
*/
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