#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

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
parameter_gamma(const double G,	const double I0, 
	const double R, const double kf);

double
calculate_c(const double gamma, const double alpha_0, const double alpha);

int
calculate_f_tide(double f_tide[9], const double G, const double m2,
	const double tilde_x[3]);

int
calculate_g(double g[9], const double G, const double m2, 
	const double alpha_0, const double alpha, const double tilde_x[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool tidal);

int
calculate_f_cent(double f_cent[9], const double omega[3]);

int
calculate_b(double b[9], const double G, const double m2, 
	const double gamma, const double alpha_0, const double alpha,
	const double tilde_x[3], const double omega[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal);

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9]);

int
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3]);

int
calculate_omega(double omega[3], const double omega_seed[3], const double G, 
	const double m2, const double I0, const double gamma, const double alpha_0, 
	const double alpha, const double tilde_x[3], const double l[3],
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal);

int
total_angular_momentum(double l_total[3],
	const double m1, const double m2,
	const double tilde_x[3], const double tilde_x_dot[3],
	const double l[3]);

double
calculate_J2(const double m, const double R, const double I[9]);

double
calculate_C22(const double m, const double R, const double I[9]);

#endif
