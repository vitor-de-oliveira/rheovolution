#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>

#include <gsl/gsl_errno.h>

#include "algelin3d.h"

/* field for 1 extended body and 1 point mass */
int
field_1EB1PM(double t, const double y[], double f[],
       		 void *params);

int
hat_map(double x_hat[9], const double x[3]);

int
construct_traceless_symmetric_matrix(double M[9], 
	const double M_main_elements[5]);

int
get_main_elements_traceless_symmetric_matrix(double M_main_elements[5], 
	const double M[9]);

double
parameter_gamma_homogeneous_body(const double G,
	const double I0, const double R);

double
calculate_c(const double gamma, const double alpha_0, const double alpha);

int
calculate_f_tide(double f_tide[9], const double G, const double m2,
	const double tilde_x[3]);

int
calculate_g(double g[9], const double G, const double m2, 
	const double alpha_0, const double alpha, const double tilde_x[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[]);

int
calculate_b(double b[9], const double G, const double m2, 
	const double gamma, const double alpha_0, const double alpha,
	const double tilde_x[3], const double omega[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[]);
int
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3]);

int
calculate_omega(double omega[3]);

#endif
