#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/* field for 1 extended body and 1 point mass */
int
field_1EB1PM(double t, const double y[], double f[],
       		 void *params);

int
hat_map(double x_hat[9], const double x[3]);

double
parameter_gamma();

int
calculate_f  (double f[9], double omega_hat[9],
              double x[3], double G, double m2);

int
calculate_lambda(double lambda[9], double f[9],
                 double b0_matrix[9], double u[9],
                 double gamma, double alpha_0,
                 double alpha);

int
calculate_b (double b[9], double f[9], double lambda[9],
             double b0_matrix[9],
             double gamma, double alpha_0);

int
calculate_l_hat (double l_hat[9], double omega_hat[9], double b[9],
		         double I0);

int
calculate_omega_hat (double omega_hat[9], double b[9], double l_hat[9],
	                 double I0);
      
#endif
