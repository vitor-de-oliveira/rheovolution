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

int
calculate_b(double b[9]);

int
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3]);

int
calculate_omega();

#endif
