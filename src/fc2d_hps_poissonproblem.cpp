#include "fc2d_hps_poissonproblem.hpp"

fc2d_hps_poisson_problem::fc2d_hps_poisson_problem(int problem_ID, double x_lower, double x_upper, double y_lower, double y_upper) :
	ID(problem_ID), x_lower(x_lower), x_upper(x_upper), y_lower(y_lower), y_upper(y_upper)
		{}

double fc2d_hps_poisson_problem::u(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 1.0;
		case LINEAR:
			return x + y;
		case LAPLACE:
			return y*sin(2*M_PI*x) + x*cos(2*M_PI*y) + 4;
		case QUAD:
			return x*x + y*y + 2*x*y;
		case POLY:
			return y*y*x*x*x*x;
		case TRIG:
			return sin(2*M_PI*x) * sin(2*M_PI*y);
		default:
			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
	}
}

double fc2d_hps_poisson_problem::f(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 0.0;
		case LAPLACE:
			return 0.0;
		case QUAD:
			return 4.0;
		case POLY:
			return 2.0*x*x*(6.0*y*y + x*x);
		case TRIG:
			return -2.0*pow(2.0*M_PI,2)*u(x,y);
		default:
			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
	}

}

double fc2d_hps_poisson_problem::dudx(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 1.0;
		case QUAD:
			return 2.0*x + 2.0*y;
		case POLY:
			return 4.0*y*y*x*x*x;
		case TRIG:
			return 2.0*M_PI*cos(2.0*M_PI*x) * sin(2.0*M_PI*y);
		default:
			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
	}
}

double fc2d_hps_poisson_problem::dudy(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 1.0;
		case QUAD:
			return 2.0*y + 2*x;
		case POLY:
			return 2.0*y*x*x*x*x;
		case TRIG:
			return sin(2*M_PI*x) * 2.0*M_PI*cos(2*M_PI*y);
		default:
			throw std::invalid_argument("[fc2d_hps_poisson_problem::u] Invalid problem_ID");
	}
}